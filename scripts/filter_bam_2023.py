from pathlib import Path
import pandas as pd
import click
import pysam
from typing import List, Dict, Tuple, Set, Union, Literal, Any
from loguru import logger
from datetime import datetime
from tqdm import tqdm
import re

__version__ = '0.0.1'

DATETIME = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
LOGFILE = Path(__file__).name.replace(".py", f"{DATETIME}.log")

IS_PRIMARY = ~0x100  ## remove secondary flag
IS_SECOUNDARY = 0x100
IS_SUPPLEMENTARY = 0x800


# # BAM= 'data/20230712_refactored/mapping/FAU74032.E6E7_66_2.refactored_ref.mmi.bam'
# BAM = Path('data/20230712_refactored/mapping/FAU74032.E6E7_66_2.refactored_ref.mmi.name_sorted.bam')
# # BAM = Path('data/20230712_refactored/mapping/FAU74032.E6E7_66_2.refactored_ref.mmi.filtered.bam')

# NCRF = 'data/20230712_refactored/NoiseCancellingRepeatFinder/FAU74032.E6E7_66_2.summary.agg.csv' 
# bed_region = 'data/refgenome/automatic/refactored_ref.tsv'
# OUTBAM = Path('data/20230712_refactored/mapping/FAU74032.E6E7_66_2.refactored_ref.mmi.filtered.bam')


## TODO:
##  1. check if read span telomere region and subtelomere region - technically 
#       checking for both subtelomere and telomere region is the same as taking the min(start) and max(end) of the coords 
#       a, check if the alignment span the telomere region
#       b, check if the alignment span the subtelomere region
##  2. check if read start with correct telomeric repeat CCTAAA for P arm in forward strand and TTAGGG for Q arm in reverse strand
##  3. if 1 and 2 are true, then check if for the secondary alignment for 1 and 2
##  4. Failing 3 check supplementary alignment for 1 and 2
##  5. Minimum mapping quality???


# 00ac554c-1dc5-485e-abaf-55204ef2d238/EcoRV/0

class Initialize_NCRF:
    """Read NCRF file"""
    df = None

    def __init__(self, ncrf_file: Path):
        self.df = pd.read_csv(ncrf_file, index_col='read_id')
        logger.info(f'Loading ncrf file: {ncrf_file}')
        self.df['corr_orient'] = self.df[['has_subtelomere',
                                          'started_correctly', 'is_terminal']].all(1) #& (self.df.segments == 1)
        self.total_reads = self.df.shape[0]
        # self.df = self.df[self.df.corr_orient == True]
        self.passed_readids = self.df.index.tolist()
        logger.info(
            f'Number of reads in ncrf file: {len(self.passed_readids)}/{self.total_reads}')
        Initialize_NCRF.df = self.df  # save to class variable

class Initialize_BED:
    """Read BED file"""
    df = None
    def __init__(self, bed_file: Path):
        self.df = pd.read_csv(bed_file, sep='\t', index_col='contig')
        logger.info(f'Loading tsv file: {bed_file}')
        Initialize_BED.df = self.df  # save to class variable

def spanning_target(targets, reads):
    t1, t2 = targets
    r1, r2 = reads
    if t1 <= r1 <= t2 <= r2:
        return True  # 3'overhange
    elif r1 <= t1 <= r2 <= t2:
        return True  # 5' overhange
    elif t1 <= r1 <= r2 <= t2:
        return True  # inside
    elif r1 <= t1 <= t2 <= r2:
        return True  # 'spanning_across')
    else:
        return False

def read_orientation(rec) -> bool:
    """Check if read is in the correct orientation by ncrf read-ids"""
    return rec.query_name in Initialize_NCRF.df.index.unique()

class AlignmentsRecord:
    def __init__(self, rec: pysam.libcalignedsegment.AlignedSegment, search_region: int = 200):
        self.rec = rec
        self.query_name = rec.query_name
        self.is_unmapped = rec.is_unmapped
        self.is_primary  = all([rec.is_supplementary == False, rec.is_secondary == False, rec.is_unmapped == False]) ## check if primary alignment
        self.reference_name = rec.reference_name
        self.reference_start = rec.reference_start
        self.reference_end = rec.reference_end
        self.target: pd.Series = Initialize_BED.df.loc[self.reference_name]  # col: telomere_start, telomere_end, subtelomere_start, subtelomere_end
        self.strand = "+" if rec.is_reverse else "-"
        self.arm = self.target.arm
        self.match_pattern = None
        self.search_region = search_region
        self.passed_filter = all(self.filter_records())
        self.mapq = rec.mapping_quality
        self.read_coverage = rec.query_alignment_length/len(rec.query_sequence)
        # self.is_unique = False
        # self.is_seen = None
        self.use_read = None

    def use_this_read(self, mapq: int = None) -> None:
        """Set read as primary and set mapping quality"""
        self.use_read = True
        self.set_tags()
        if not self.is_primary:
            self.set_as_primary()
            self.set_mapping_quality(mapq)
        
    def span_telomere_region(self) -> bool:
        """Check if alignment span the telomere region"""
        return spanning_target((self.target.telomere_start, self.target.telomere_end), (self.reference_start, self.reference_end))

    def span_subtelomere_region(self) -> bool:
        """Check if alignment span the subtelomere region"""
        return spanning_target((self.target.subtelomere_start, self.target.subtelomere_end), (self.reference_start, self.reference_end))

    def check_alignment_orientation(self) -> bool:
        """Check if alignment is in the correct orientation"""
        # TODO: check if the sequence alway in the foward orientation by need to reverse complement!!
        # Read orientation from pysam: "The sequence is returned as it is stored in the BAM file.
        # (This will be the reverse complement of the original read sequence if the mapper has aligned the read to the reverse strand.)"
        # https://pysam.readthedocs.io/en/latest/api.html?highlight=query_sequence%20#pysam.AlignedSegment.query_sequence
        # Ideally, we want to check the aligned sequence but in the cases with telomeres are masked, we can't do that
        ## TODO: FIXES THIS PART OF THE CODE-- DONE
        if self.arm == 'P':  # left
            pattern = "CCCTAA"
        elif self.arm == "Q":  # right
            pattern = "TTAGGG"
        elif self.arm == "M":
            return False
        else:
            raise ValueError(f"Arm must be either P, Q oe M")
        self.match_pattern = re.compile(".*?"+pattern)
        res = self.match_telomeric_repeat()
        # logger.debug(
            # f'{self.query_name}: Arm: {self.arm} Strand: {self.strand}\nPattern:{pattern} MatchEnd {"Start" if self.arm =="P" else "End"} - Match: {res}\nStart: {self.rec.query_sequence[:200]}\nEnd: {self.rec.query_sequence[-200:]}')
        return True if res else False

    def match_telomeric_repeat(self) -> bool:
        """Check if the alignment start with the correct telomeric repeat"""
        if self.arm == "P": # left or start of seq
            return self.match_pattern.match(self.rec.query_sequence[:self.search_region])
        elif self.arm == "Q": # right or end of seq
            return self.match_pattern.match(self.rec.query_sequence[-self.search_region:])
        else:
            raise ValueError(f"Arm must be either P or Q")

    def filter_records(self) -> Tuple[bool, bool, bool]:
        """Filter records"""
        return self.span_telomere_region(), self.span_subtelomere_region(), self.check_alignment_orientation()
    
    @property
    def atype(self):
        """ALignment type"""
        if self.is_primary:
            return "primary"
        elif self.rec.is_secondary:
            return "secondary"
        elif self.rec.is_supplementary:
            return "supplementary"
        elif self.rec.is_unmapped:
            return "unmapped"
        else:
            # logger.warning(f"Unknown alignment type: {self.query_name} {self.atype}  (primary, sec, sup, unmap) ({self.primary}, {self.rec.is_secondary:}, {self.rec.is_supplementary}, {self.rec.is_unmapped})")
            return "unknown"

    def set_mapping_quality(self, mapq: int):
        self.rec.mapping_quality = mapq

    def set_as_primary(self):
        """Set alignment as primary"""
        self.rec.flag = self.rec.flag & IS_PRIMARY
        logger.info(f"Set {self.query_name} {self.reference_name} {self.reference_start} {self.reference_end} as primary alignment")

    def set_as_secondary(self):
        """Set alignment as secondary"""
        self.rec.flag = self.rec.flag | IS_SECOUNDARY
        self.rec.mapping_quality = 0
        logger.info(f"Set {self.query_name} {self.reference_name} {self.reference_start} {self.reference_end} as primary alignment")

    def set_tags(self):
        self.rec.set_tag('AM', self.arm)  ## set arm
        self.rec.set_tag("CM", self.target.chromosome, replace=True)  ## set chromosome
        self.rec.set_tag('AP', 1, replace=True)  ## passed filter 

    def header(self):
        return "\t".join(["query_name", "atype", "seq_length", 
                          'flag', 'mapq', 'read_coverage',
                          "reference_name", "reference_start", "reference_end", 
                          "strand", "span_telomere_region", "span_subtelomere_region", 
                          "check_alignment_orientation", "use_read"]
                          )+ "\n"

    def to_tsv(self):
        return "\t".join(map(str, 
                             [self.query_name, self.atype, len(self.rec.query_sequence), 
                              self.rec.flag, self.mapq, self.read_coverage,
                              self.reference_name, self.reference_start, self.reference_end, 
                              self.strand, self.span_telomere_region(), self.span_subtelomere_region(),
                              self.check_alignment_orientation(), self.use_read])
        )+'\n'
    

atype_order = {
    "primary": 0,
    "secondary": 1,
    "supplementary": 2,
    "unmapped": 3
}


def set_option_for_passed_reads(rec):
    ...



def read_bam(bamfile_handler) -> AlignmentsRecord:
    """Read **name sorted** bam files and yield the reads that passed the filter.
    Returns the reads in the following order: primary, secondary, supplementary as AlignmentRecord objects.
    Th
    
    """
    # readid_passed_tracker = []
    # primary_is_seen = {} ## keep track of primary alignment
    read_has_passed = {} ## keep track of reads that has passed filter
    reads = []
    current_read = None 

    missing_reads = ["4a594caf-b2e2-48d2-a8d2-d5bde90f6f8d/EcoRV/0", "4a2c5138-457c-4c60-824b-824c951e75b6/EcoRV/0",
                     "49ca673f-dccb-480c-8c5d-8a5af355d36a/EcoRV/0"]


    for record in bamfile_handler:
        if record.is_unmapped:
            continue  # skip unmapped reads
        if read_orientation(record) == False: ## skip reads that are not in the ncrf file
            if record.query_name in missing_reads:
                raise ValueError(f"Read orientation is False: {record.query_name} {record.reference_name} \n{record.rec.query_sequence}")
            continue # skip reads that are not in the ncrf file
        if not record.reference_name in Initialize_BED.df.index.unique(): 
            continue
        r = AlignmentsRecord(record)
        # logger.debug(f'{r.query_name} {r.atype}: {r.is_primary} {r.passed_filter} {r.reference_name}')
        
        if current_read is None:
            current_read = r.query_name
            reads.append(r)
            continue
        elif current_read == r.query_name:
            reads.append(r)
            continue
        elif current_read != r.query_name:
            ## sort the reads by alignment type and yield the reads
            reads.sort(key=lambda x: atype_order[x.atype])
            # max_read_coverage = max([x.rec.query_alignment_length/len(x.rec.query_sequence) for x in reads])
            logger.debug(f'Sorting by alignment types: {[x.atype for x in reads]}')

            found_passed_read = False
            mapq = reads[0].rec.mapping_quality
            for rec in reads:
                ## read is already passed filter 
                if found_passed_read or read_has_passed.get(rec.query_name, False):
                # if rec.query_name in readid_passed_tracker or rec.query_name in primary_is_seen:
                    logger.debug(f'{rec.query_name} {rec.atype}: Already passed filter')
                    rec.use_read = False
                    yield rec
                if not found_passed_read:    
                    if rec.is_primary and rec.passed_filter:
                        logger.info(
                            f'{rec.query_name}: (telo, subtelo, orient, atype, ref) {rec.filter_records()} {rec.atype} {rec.reference_name}')
                        # readid_passed_tracker.append(rec.query_name)
                        read_has_passed[rec.query_name] = True
                        # rec.use_read = True
                        rec.use_this_read()
                        found_passed_read = True

                    elif rec.is_primary and not rec.passed_filter:
                        rec.set_as_secondary() ## set as secondary 

                        rec.use_read = False

                    elif not read_has_passed.get(rec.query_name, False) and rec.rec.is_secondary and rec.passed_filter:
                        logger.info(
                            f'{rec.query_name}: (telo, subtelo, orient, atype, ref) {rec.filter_records()} {rec.atype} {rec.reference_name}')
                        # readid_passed_tracker.append(rec.query_name)
                        read_has_passed[rec.query_name] = True
                        ## set mapping quality to the same as primary
                        rec.use_this_read(mapq)
                        # rec.use_read = True
                        # rec.set_as_primary()
                        # rec.set_mapping_quality(mapq) 
                        found_passed_read = True

                    elif not read_has_passed.get(rec.query_name, False) and rec.rec.is_supplementary and rec.passed_filter:
                        logger.info(
                            f'{rec.query_name}: (telo, subtelo, orient, atype, ref) {rec.filter_records()} {rec.atype} {rec.reference_name}')
                        # readid_passed_tracker.append(rec.query_name)
                        read_has_passed[rec.query_name] = True
                        rec.use_read = True
                        found_passed_read = True
                    else:
                        rec.use_read = False ## don't use the read, failed the filter
                yield rec
            
            reads = [r]  ## add the current read to the list
            current_read = r.query_name  ## update the current read
 

def pass_filters(rec: AlignmentsRecord) -> bool:
    """Check if read has already passed filter"""
    if all([rec.use_read == True, rec.passed_filter == True]): ## explicit check
        return True
    return False


@click.command()
@click.argument("bamfile")
@click.argument("ncrf")
@click.argument("tsv")
@click.argument("bamout")
@click.option("--debug", '-d', is_flag=True, help="Debug mode")
@click.option("--skip-supplementary", '-s', is_flag=True, help="Skip supplementary alignment from output bam")
@click.version_option(__version__)
def main(bamfile: str, ncrf: str, tsv: str, bamout: str, debug: bool, skip_supplementary: bool):
    """Main function"""
    bamout = Path(bamout)
    if not bamout.parent.exists():
        bamout.parent.mkdir(parents=True)
    LOGFILE = bamout.with_suffix(f".{DATETIME}.log")

    logger.add(LOGFILE)
    logger.info(f"Logfile: {LOGFILE}")
    logger.info(f'Script version: {__version__}')
    # if debug:
    # logger.add(LOGFILE, level="DEBUG")
    logger.info("Debug mode")
    logger.info(f"Bamfile: {bamfile}")
    logger.info(f"NCRF: {ncrf}")
    logger.info(f"TSV: {tsv}")
    logger.info(f"BAMOUT: {bamout}")
# logger.
    ncrf = Initialize_NCRF(ncrf)
    anchor_bed = Initialize_BED(tsv)


    # missing_reads = ["4a594caf-b2e2-48d2-a8d2-d5bde90f6f8d/EcoRV/0", "4a2c5138-457c-4c60-824b-824c951e75b6/EcoRV/0",
    #  "49ca673f-dccb-480c-8c5d-8a5af355d36a/EcoRV/0"]


    # for i in missing_reads:
    #     # raise ValueError(ncrf.df.index.to_list())
    #     if i in ncrf.df.index.to_list():
    #         print(f'>>> {i} YES')
    #     else:
    #         print(f'{i} Not found in ncrf file')
    # raise ValueError(
    #     ncrf.df.loc[['4a594caf-b2e2-48d2-a8d2-d5bde90f6f8d/EcoRV/0']])

    total_passed = 0
    header_writen = False
    bam_fh = pysam.AlignmentFile(bamfile, 'rb')
    out = pysam.AlignmentFile(bamout, 'wb', template=bam_fh)

    # with pysam.AlignmentFile(BAM, 'rb') as bam_fh, pysam.AlignmentFile(OUTBAM, 'wb', template=bam_fh) as out, open(OUTBAM.with_suffix('.records.tsv'), 'w') as tsv_out:
    with open(bamout.with_suffix('.tsv'), 'w') as tsv_out:
        for rec in read_bam(bam_fh):
            # if total_passed == 100:
            #     break
            # if rec.query_name in missing_reads:
                # raise ValueError(f"Missing reads {rec.query_name} {rec.filter_records()} {rec.is_primary} {rec.passed_filter} {rec.reference_name} \n{rec.rec.query_sequence}")
            if not header_writen:
                tsv_out.write(rec.header())
                header_writen = True
            tsv_out.write(rec.to_tsv())
            # if rec.query_name.startswith("8647905c-77f8-4f12-a354-0e355d8a1a1b"):
            #     raise ValueError(f"{rec.filter_records()} {rec.is_primary} {rec.passed_filter} {rec.reference_name} \n{rec.rec.query_sequence}")
            if pass_filters(rec):
                if skip_supplementary and rec.rec.is_supplementary:
                    continue
                out.write(rec.rec)
                logger.info(f"Read passed: {rec.query_name} {rec.atype} {rec.reference_name} {rec.passed_filter} {rec.use_read}")
                total_passed += 1
    bam_fh.close()
    out.close()
            # print(rec.query_name, reference_name)
    logger.info(f"Number of reads passed: {total_passed}")
    logger.info(f"Filtered bam file written to: {bamout}")
        # break

if __name__ == '__main__':
    main()
