
from pathlib import Path
from dataclasses import dataclass
import click
import regex
from intervaltree import Interval, IntervalTree
import pysam
from collections import Counter 
from typing import List, Dict, Tuple, Set, Union, Literal, Any 
import wget 
from loguru import logger
from datetime import datetime
from tqdm import tqdm
import re

__version__ = '0.0.1'

DATETIME = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
LOGFILE =  Path(__file__).name.replace(".py", f"{DATETIME}.log")
logger.add(LOGFILE)
logger.info(f'"Logfile: {LOGFILE}")')
logger.info(f'Script version: {__version__}")')


def telomere_pattern() -> regex.Regex:
    """regex telomere pattern, not case sensitive, 1 error (ins, del, sub) allowed """
    pattern = '(TTAGGG|CCCTAA){e<=%s}' % 0  # allow for 1 error (ins, del, sub)
    return regex.compile(pattern, flags=regex.IGNORECASE)

def get_telomere_region(seq: str) -> IntervalTree:
    """Get the telomere region"""
    tree = IntervalTree()
    matched_motif = []
    data = []
    # Find telomere pattern matches in the sequence
    t_pattern = telomere_pattern()
    for res in t_pattern.finditer(seq, concurrent=True, timeout=1, overlapped=True):
        data.append(res)
        matched_motif.append(res.group())
    # Count matched motif
    matched_motif_counts = Counter(matched_motif)
    common_kmer = matched_motif_counts.most_common(1)[0][0]
    # Add intervals to interval tree based on the common kmer
    for res in data:
        if res.group() == common_kmer:
            tree.addi(res.start(), res.end(), common_kmer)
    return tree

def merge_tree(tree: IntervalTree, max_break_point: int = 250) -> IntervalTree:
    """Merge overlapping intervals with X bp"""
    # Merge the intervals if the distance is less than 20
    def reducer(old, new):
        if old is None or new is None:
            return '.'
        if old == new:
            return new
        return f"{old}, {new}"

    merged_tree = IntervalTree()
    data = set()
    for interval in tree:
        start = interval.begin
        end = interval.end
        # Find overlapping intervals within a distance 
        overlapping_intervals = tree.overlap(start - max_break_point, end + max_break_point)
        # Calculate the new merged interval
        data = set([interval.data for interval in overlapping_intervals])
        if len(data) == 1:
            data = list(data)[0]
        else:
            raise Exception(f"Data missing: {data}")
        merged_start = min([interval.begin for interval in overlapping_intervals])
        merged_end = max([interval.end for interval in overlapping_intervals])
        # Create the merged interval
        merged_interval = Interval(merged_start, merged_end, data=data)
        # Add the merged interval to the new IntervalTree
        merged_tree.add(merged_interval)

    merged_tree.merge_overlaps(strict=False, data_reducer=reducer)
    return merged_tree

def adjust_coords(reference_length: int, search_region: int, tree: IntervalTree) -> IntervalTree:
    """Adjust the coordinates to the reference genome"""
    def _adjust_coord(reference_length: int, search_region: int, coord) -> int:
        """adjusted_coord = (reference_length - search_region) + coord"""
        return (reference_length - search_region) + coord
    new_tree = IntervalTree()
    for interval in tree:
        new_tree.addi(
            _adjust_coord(reference_length, search_region, interval.begin), 
            _adjust_coord(reference_length, search_region, interval.end),
            data=interval.data
            )
    return new_tree

def get_telomere_length(tree: IntervalTree) -> int:
    """Get the telomere length"""
    return sum([interval.end - interval.begin for interval in tree])

def filter_for_arm(tree: IntervalTree, is_q: bool=False) -> IntervalTree:
    """Filter for arm"""
    sorted_tree = sorted(tree)
    if is_q:
        return list(sorted_tree)[-1]
    else:
        # is p arm
        return list(sorted_tree)[0]
    
def subtelomere_coords(tree: IntervalTree, search_region: int, is_q: bool = False) -> IntervalTree:
    """Get the subtelomere coordinates"""
    if is_q:
        # is q arm on ajdusted coordinates
        # return Interval((tree.begin-search_region)+tree.length(), tree.begin)
        return Interval(0, tree.begin)
    else:
        # is p arm
        return Interval(tree.end, search_region)

def mask_region(seq: str, start: int, end: int) -> str:
    """Mask the region with Ns"""
    return seq[:start] + 'N'*(end-start) + seq[end:]

def get_chrom_name(config: str) -> str:
    """Get the chromosome name"""
    r = re.search(r'(chr[0-9XY]+)', str(config), flags=re.IGNORECASE)
    if r:
        return r.group(1)
    return None

class Genome:
    def __init__(self, url: str, filename: str, genome, keep_middle: bool = False, output_dir: Path = None, skip_trimming: bool = False) -> None:
        self.url = url
        self.keep_middle = keep_middle
        self.output_dir = self._output_dir(output_dir)
        self.filename = self.output_dir / filename
        self.genome = genome  ## genome name i.e. HG002 or CHM13
        self.skip_trimming = skip_trimming

    def _output_dir(self, output_dir: Path) -> None:
        """Create the output directory if it doesn't exist"""
        if output_dir is None:
            output_dir = Path('refgenome')
        if not output_dir.exists():
            output_dir.mkdir(parents=True, exist_ok=True)
        return output_dir

    def download_ref(self) -> str:
        """Download the reference genome"""
        if not self.filename.exists():
            logger.info(
                f'Please be patient, downloading {self.url} to {self.filename}')
            self.filename = wget.download(self.url, out=str(self.filename))
            logger.info(f'Downloaded: {self.filename}')
        else:
            logger.info(f'{self.filename} already exists, skipping download ...')
        return self.filename

@dataclass
class Record:
    chrom: str ## chromsome name
    config: str  ## orgional config
    arm: Literal['P', 'Q', "M"]
    seq: str
    genome: Genome
    telomere_coords: IntervalTree = None
    subtelomere_coords: IntervalTree = None
    seq_coords: IntervalTree = None

    def __post_init__(self):
        self.arm = self.arm.upper()
        if self.arm not in ['P', 'Q', 'M']:
            raise ValueError(f"Invalid arm: {self.arm}")

    @property
    def header(self) -> str:
        """Return the header 
        FORMAT: chr.arm|genome|config:start-end
        """
        if self.arm == 'M':
            return f'{self.chrom}.{self.arm}|{self.genome.genome}|{self.config}:None-None'
        elif self.arm == "Q":
            return f'{self.chrom}.{self.arm}|{self.genome.genome}|{self.config}:{self.seq_coords.end}-{self.seq_coords.begin}'
        return f'{self.chrom}.{self.arm}|{self.genome.genome}|{self.config}:{self.seq_coords.begin}-{self.seq_coords.end}'

    def format_sequence(self) -> str:
        """Fasta formated sequence"""
        return f'>{self.header}\n{self.seq}\n'

    def to_tsv(self) -> str:
        """Write the telomere/sub coordinates to a tsv file
        FORMAT: header, telomere_start, telomere_end, telomere_length, telomere_common_kmer, subtelomere_start, subtelomere_end, chromosome, arm
        """
        values = [self.header, self.telomere_coords.begin, self.telomere_coords.end, self.telomere_coords.length(), self.telomere_coords.data.upper(), 
                   self.subtelomere_coords.begin, self.subtelomere_coords.end, self.chrom, self.arm]
        if self.arm == 'M':
            values = [self.header, 'None', 'None', 'None', 'None', 'None', 'None', self.chrom, self.arm]
        return '\t'.join(list(map(str,values))) + '\n'

    @staticmethod
    def tsv_header():
        header = ['contig', 'telomere_start', 'telomere_end', 'telomere_length', 'telomere_common_kmer']
        header += ['subtelomere_start', 'subtelomere_end', 'chromosome', 'arm']
        return '\t'.join(list(map(str,header))) + '\n'


def stong_chrom_arm_name(config: str) -> Tuple[str, str]:
    r = re.search(r'(^[0-9XY]+)(P|Q)', str(config), flags=re.IGNORECASE)
    if r:
        return f"chr{r.group(1)}", r.group(2)
    return None

def process_stong_etal_assembly(rec: Any, genome: Genome, search_region: int) -> Record:
    """Process the Strong et al. assembly. 
    Note: The sequences of stong et al assemblies are oriented with the telomere on the left and aligned
    Ref: Stong, N., et al 2014. 
    Subtelomeric CTCF and cohesin binding site organization using improved subtelomere assemblies and a novel annotation pipeline. 
    Genome Research"""
    chrom, arm = stong_chrom_arm_name(rec.name)
    if chrom is None:
        raise Exception(f"Unable to parse the chromosome name: {rec.name}")
    if arm is None:
        raise Exception(f"Unable to parse the arm name: {rec.name}")
    ## all sequences are oriented with the telomere on the forward (left)
    seq = rec.sequence[:search_region]
    return Record(
            chrom=chrom,
            config=rec.name,    
            arm=arm.upper(),
            seq=rec.sequence[:search_region],
            genome=genome,
            telomere_coords=Interval(1, 50, data= "TTAGGG" if arm.upper() == "Q" else 'CCCTAA'), ## TODO: get the telomere coordinates
            subtelomere_coords=Interval(50, len(seq)),
            seq_coords=Interval(0, len(seq))
            )


def process_genome(genome: Genome, search_region: int = 25_000, max_break_point: int = 200, hg002_pat: bool = False) -> None:
    """Process the genome"""
    # ref = Path('/mmfs1/data/active/projects/202009_salk_telomere/ANALYSIS/prughani/refgenome/HG002/HG002_assembly.v0.7.fasta')
    # assert genome.filename.exists()

    if not genome.filename.exists():
        genome.download_ref()
        
    logger.info(f"Only keeping MATERNAL genome for HG002: {hg002_pat}")
    logger.info(f"Processing the genome {genome.filename}")
    
    # search_region = search_region
    # max_break_point = max_break_point ## merge overlapping intervals within 200 bp distance
    data = []
    for rec in tqdm(pysam.FastxFile(genome.filename), desc=f"Searching for telomeres: {genome.filename}"):
        # if not rec.name.startswith('chr') or rec.name in ["chrEBV", 'chrM']:

        if genome.genome == "HG002" and hg002_pat:
            if not re.match('(chr[0-9XY]{1,2}_PATERNAL)', rec.name, flags=re.IGNORECASE):
                logger.debug(f"Skipping contig: {rec.name}")
                continue

        ## process the strong et al. assembly differently, include the whole chromosome as the subtelomere region
        if genome.genome == 'strong_etal_2014' or genome.skip_trimming is True:
            data.append(process_stong_etal_assembly(rec, genome, search_region))
            continue

        if get_chrom_name(rec.name) is None:
            continue

     
        parm = rec.sequence[:search_region]
        qarm = rec.sequence[-search_region:]
    
        ptree = get_telomere_region(parm)
        qtree = get_telomere_region(qarm)
        ## Merge the overlapping intervals within `max_break_point` bp distance
        mptree = merge_tree(ptree, max_break_point=max_break_point)
        mqtree = merge_tree(qtree, max_break_point=max_break_point)

        ## Adjust the coordinates to the reference genome for q arm
        adj_qtseq_coord = filter_for_arm(adjust_coords(len(rec.sequence), len(qarm), mqtree), is_q=True)
        adj_qtsubt_coord = subtelomere_coords(adj_qtseq_coord, len(qarm), is_q=True)

        pcoord = filter_for_arm(mptree, is_q=False)
        qcoord = filter_for_arm(mqtree, is_q=True)

        subtpcoord = subtelomere_coords(pcoord, len(parm), is_q=False)
        subtqcoord = subtelomere_coords(qcoord, len(qarm), is_q=True)
        
        data.append(Record(
            chrom=get_chrom_name(rec.name),
            config=rec.name,    
            arm='P',
            seq=parm,
            genome=genome,
            telomere_coords=pcoord,   
            subtelomere_coords=subtpcoord,
            seq_coords=Interval(0, len(parm))
            )
        )
        
        data.append(Record(
            chrom=get_chrom_name(rec.name),
            config=rec.name,    
            arm='Q',
            seq=qarm,
            genome=genome,
            telomere_coords=qcoord,   
            subtelomere_coords=subtqcoord,
            seq_coords=Interval(len(rec.sequence)-len(qarm), len(rec.sequence))

            )
        )
        # raise ValueError(f"{data[-1].seq_coords}" )
        if genome.keep_middle:
            # Keep the middle region
            middle = rec.sequence[search_region:-search_region]
            data.append(Record(
                chrom=get_chrom_name(rec.name),
                config=rec.name,    
                arm='M',
                seq=middle,
                genome=genome,
                telomere_coords=None,   
                subtelomere_coords=None,
                seq_coords=Interval(len(parm), len(rec.sequence)-len(qarm))
                )
            )
    
    return data


@click.command()
@click.argument('prefix', type=str)
@click.option('--search_region', '-s', default=25_000, type=int, help='Search region for telomere')
@click.option('--max_break_point', '-x', default=200, type=int, help='Max break point for merging overlapping intervals')
@click.option('--output_dir', '-o', default='refgenome', type=click.Path(file_okay = False), help='Output directory')
@click.option('--working_dir', '-w', default='workdir', type=click.Path(file_okay = False), help='Working directory')
@click.option('--hg002_pat', '-p', is_flag=True, help='Output the paternal genome for HG002')
# @click.option('--coverant','-c', default=False, is_flag=True, help='Convert to full legnth genomic coordinates')
@click.version_option(version=__version__)
def main(prefix, output_dir, search_region, max_break_point, working_dir, hg002_pat):
    """
    Generate a reference genome with telomere and subtelomere regions for HG002, CHM13, and Strong et al. assemblies. 
    Only middle region for CHM13 assembly is kept.
    \n
    Params:
        Prefix: output file
    \n
    Outputs:
        1. fasta file with the telomere and subtelomere seqeunces
        2. tsv file with the telomere and subtelomere coordinates without middle genome region for HG002 and stong et al. assemblies
    \n
    Example:
    python bam_filter_telomere_v2.py --search_region 25_000 --max_break_point 200 --output_dir refgenome --working_dir workdir my_refgenome

    """
    working_dir = Path(working_dir)

    hg002 = Genome(
        url='https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/assemblies/drafts/assembly.v0.7.fasta',
        filename='HG002_assembly.v0.7.fasta',
        genome='HG002',
        keep_middle=False,
        output_dir=working_dir
        )
        
    chm13 = Genome(
        url="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz",
        filename='chm13v2.0.fa.gz',
        genome='CHM13',
        keep_middle=True,
        output_dir=working_dir
        )

    strong = Genome(
        url='https://genome.cshlp.org/content/suppl/2014/04/16/gr.166983.113.DC1/Supplemental_FileS1.txt',
        filename='strong_etal_2014.subtelomere.fasta',
        genome='strong_etal_2014',
        keep_middle=False,
        output_dir=working_dir,
        skip_trimming=True,
        )

    genomes = [hg002, chm13]
    if hg002_pat:
        # hg002.skip_trimming = True
        hg002.keep_middle = True
        genomes = [hg002]
    output = Path(output_dir)
    if not output.exists():
        output.mkdir(parents=True, exist_ok=True)
    output_file = output / f"{prefix}.fasta"
    with open(output_file, 'w') as fasta, open(output_file.with_suffix('.tsv'), 'w') as tsv:
        tsv.write(Record.tsv_header())
        for genome in genomes:
            for record in tqdm(process_genome(genome, search_region=search_region, max_break_point=max_break_point, hg002_pat=hg002_pat)):
                fasta.write(record.format_sequence())
                if record.arm != 'M': ## only write the P and Q arm coordinates
                    tsv.write(record.to_tsv())
    logger.info(f"Done! Output file: {output_file}")

if __name__ == "__main__":
    main()
