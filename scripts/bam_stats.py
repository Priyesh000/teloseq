#!/usr/bin/env python
import pandas as pd
import pysam
from pathlib import Path
from collections import OrderedDict

# from multiprocessing import Pool
from functools import partial
import concurrent.futures
import sys
import logging
import argparse
import re


def get_args():
    parser = argparse.ArgumentParser("Get stats from bam AlignmentFile")
    parser.add_argument("bam", metavar="FILE", help="Path to bam")
    parser.add_argument(
        "--threads",
        "-t",
        metavar="INT",
        default=4,
        type=int,
        help="number of threads. Default[4]",
    )
    parser.add_argument(
        "--tol",
        "-x",
        metavar="INT",
        default=50_000,
        type=int,
        help="tolerance for chromosome end classification. Default[50,000]",
    )
    parser.add_argument("--bed", "-b", metavar="FILE", default=None, help="Bed FILE")
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="increase output verbosity"
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--stdout",
        "-s",
        action="store_true",
        help="Write to standard out (STDOUT). Default[False]",
    )
    group.add_argument("--output", "-o", metavar="FILE", help="Output csv")

    parser.add_argument(
        "--cigar", "-c", action="store_true", help="parse the cigar string and return match, insetion, deletion and cliped bases [False]"
    )
    return parser.parse_args()

def at_ends(x, end, tol=50000,):
    """
    Assign the end where chromosome
    """
    if x <= tol:
        ## 0-tol
        return 5
    elif (end - tol) <= x <= end:
        return 3
    else:
        return 0

def edit_distance_cigar(rec):
    """
    use the cigar sting to calcualate the edit_distance and MD tag
    """
    match = sum(
        [len(item) for item in re.split("[0-9^]", rec.get_tag("MD"))]
    )  # Parse MD string to get mismatches/deletions
    insertion = sum(
        [item[1] for item in rec.cigartuples if item[0] == 1]
    )  # Parse cigar to get insertions
    return (match + insertions) / rec.query_alignment_length

def edit_distance(rec):
    try:
        NM = rec.get_tag("NM")
        return NM / rec.query_alignment_length
    except:
        # return edit_distance_cigar(rec)
        return None


def parse_cigar(cigartuples):
    ## https://pysam.readthedocs.io/en/latest/api.html

    # cig_values = {

    #     0: "matches",
    #     1: "insertions",
    #     2: "deletions",
    #     3: "ref_skip",
    #     5: "soft_clip",
    #     6: "hard_clip"
    # }

    operation_values = OrderedDict([
        (0, 0),  # match
        (1, 0),  # insertions
        (2, 0),  # deletions
        (3, 0),  # ref_skip
        (4, 0),  # soft_clip
        (5, 0),   # hard_clip
        ('left_soft_clip', 0),
        ('right_soft_clip', 0),
        ('identity', 0)
    ])
    cig_len = len(cigartuples)

    for idx, (operation, length) in enumerate(cigartuples):
        if operation == 4:
            if idx > (cig_len/2):
                operation_values['right_soft_clip'] += length
            else:
                operation_values['left_soft_clip'] += length
        if operation <= 5:
            operation_values[operation] += length
    # return {cig_values[k]: v for k, v in operation_values}
    operation_values['identity'] = (
            100 * (operation_values[0] / sum([operation_values[x] for x in range(3)])
                )
            )
    return operation_values

def bam_mapping_stats(bam, region, parse_cigar_string , tol=50000):
    data = []

    bam_handle = pysam.AlignmentFile(
        bam, "rb"
    )  ### open a new bam in mem for each process
    if isinstance(region, tuple):
        region = {k: v for k, v in zip(["contig", "start", "stop"], region)}

        # handler = bam_handle.fetch(*region, multiple_iterators=True)
    else:
        # handler = bam_handle.fetch(contig=region, multiple_iterators=True)
        region = {"contig": region}
    for rec in bam_handle.fetch(**region, multiple_iterators=True):
        if rec.is_unmapped:
            continue

        if rec.is_supplementary:
            atype = "supplementary"
        elif rec.is_secondary:
            atype = "secondary"
        else:
            atype = "primary"
        tags = dict(rec.tags)


        d = [
            rec.query_name,
            rec.query_alignment_start,
            rec.query_alignment_end,
            rec.query_alignment_length,
            rec.infer_read_length(),
            "-" if rec.is_reverse else "+",
            rec.reference_name,
            rec.reference_start,
            rec.reference_end,
            bam_handle.get_reference_length(rec.reference_name),
            rec.mapping_quality,
            atype,
            at_ends(
                rec.reference_start,
                bam_handle.get_reference_length(rec.reference_name),
                tol=tol,
            ),
            rec.get_tag("NM") / rec.query_alignment_length,
            tags.get("AM",-1),
            tags.get("AP",-1),
            tags.get("RG",-1),
        ]
        if parse_cigar_string:
            d = d + list(parse_cigar(rec.cigartuples).values())

        data.append(d)
    return data

def colnames(parse_cigar_string):
    col_names =  [
        "read_id",
        "qstart",
        "qend",
        "q_align_len",
        "read_len",
        "strand",
        "rname",
        "rstart",
        "rend",
        "ref_length",
        "qmap",
        "atype",
        "strand_end",
        "edit_distance",
        "arm",
        "anchored",
        "read_group"

    ]
    if parse_cigar_string:
        cig_cols = [
            "matches",
            "insertions",
            "deletions",
            "ref_skip",
            "soft_clip",
            "hard_clip",
            'left_soft_clip',
            'right_soft_clip',
            'identity'
        ]
        return col_names + cig_cols
        
# def calcuate_accurary():
#     # Calculate global identity and accuracy:
#     base_stats['identity'] = float(
#         base_stats['match']) / (base_stats['match'] + base_stats['mismatch'])
#     base_stats['accuracy'] = 1.0 - (float(base_stats['insertion'] +
#                                           base_stats['deletion'] + base_stats['mismatch'] + base_stats['clipps']) / base_stats['aln_length'])


def get_chromosome(bam):
    return check_bam_index(bam).references

# sub = tags['NM'] - ins - delt


def process_bam(bam, bed=None, threads=2, stout=False, parse_cigar_string='',  **kwargs):

    if bed:
        regions = parse_bed(bed)
    else:
        regions = get_chromosome(bam)

    func = partial(bam_mapping_stats, bam,
                   parse_cigar_string=parse_cigar_string, ** kwargs)
    data = []
    if stout:
        sys.stdout.write(
            ",".join(colnames(parse_cigar_string=parse_cigar_string)) + "\n")
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        for region, results in zip(regions, executor.map(func, regions)):
            logging.debug(f"Processed region: {region}")
            for result in results:
                if stout:
                    sys.stdout.write(",".join(map(str, result)) + "\n")
                else:
                    data.append(result)
    if not stout:
        return pd.DataFrame(data, columns=colnames(parse_cigar_string=parse_cigar_string))
    else:
        return False

def check_bam_index(bam):
    samfile = pysam.AlignmentFile(bam, "rb")
    if not samfile.has_index():
        logging.info("Bam index is missing. Creating a new index")
        pysam.index(bam)
        return pysam.AlignmentFile(bam, "rb")  ## reload Bam
    return samfile

def file_validation(fn):
    if isinstance(fn, bool) or fn is None:
        return None
    fn = Path(fn)
    if fn.exists():
        logging.info("Output file exists and will be overwritten")
    elif not fn.parent.exists():
        fn.parent.mkdir(parents=True)
    return fn

def parse_bed(bed):
    data = []
    with open(bed) as fh:
        for line in fh:
            line = line.strip()
            attrs = line.split("\t")
            start, end = map(int, attrs[1:3])
            data.append((str(attrs[0]), start, end))
    logging.debug(data)
    return data

def main():
    args = get_args()
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO

    logging.basicConfig(
        stream=sys.stderr,
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%d/%m/%Y %H:%M:%S",
        level=level,
    )

    logging.info(args)
    if args.stdout:
        logging.info("Writting to STDOUT")
    output_fn = file_validation(args.output)
    df = process_bam(
        args.bam, bed=args.bed, stout=args.stdout, 
        threads=args.threads, tol=args.tol, parse_cigar_string=args.cigar
    )

    if isinstance(df, pd.DataFrame):
        logging.info(f"Writting to {args.output}")
        df.to_csv(output_fn, index=False)
    logging.info("Finished!")

if __name__ == "__main__":
    main()
    