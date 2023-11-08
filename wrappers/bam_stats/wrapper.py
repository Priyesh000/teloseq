from snakemake import shell


import pandas as pd
import pysam
from pathlib import Path
from collections import OrderedDict

from functools import partial
import concurrent.futures
import sys
# import logging
# import argparse
import re


snakemake = snakemake ## noqa


#'set +u; python scripts/bam_stats.py  -t {threads} --tol 50000 -s  {input} > {output}; set -u'

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


def bam_mapping_stats(bam, region, tol=50000):
    data = []
    ## open a new bam in mem for each process
    bam_handle = pysam.AlignmentFile(
        bam, "rb"
    )  
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

        d = (
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
            rec.get_tag("NM"),
            rec.get_tag("NM") / rec.query_alignment_length,
        )

        data.append(d)
    return data


def colnames():
    return (
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
        "norm_edit_distance"
    )


def get_chromosome(bam):
    return check_bam_index(bam).references


def process_bam(bam, bed=None, threads=2, stout=False, **kwargs):

    if bed:
        regions = parse_bed(bed)
    else:
        regions = get_chromosome(bam)

    func = partial(bam_mapping_stats, bam, **kwargs)
    data = []
    if stout:
        sys.stdout.write(",".join(colnames()) + "\n")
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        for region, results in zip(regions, executor.map(func, regions)):
            # logging.debug(f"Processed region: {region}")
            for result in results:
                if stout:
                    sys.stdout.write(",".join(map(str, result)) + "\n")
                else:
                    data.append(result)
    if not stout:
        return pd.DataFrame(data, columns=colnames())
    else:
        return False


def check_bam_index(bam):
    samfile = pysam.AlignmentFile(bam, "rb")
    if not samfile.has_index():
        print("Bam index is missing. Creating a new index") 
        # logging.info("Bam index is missing. Creating a new index")
        pysam.index(bam)
        return pysam.AlignmentFile(bam, "rb")  ## reload Bam
    return samfile


def file_validation(fn):
    if isinstance(fn, bool) or fn is None:
        return None
    fn = Path(fn)
    if fn.exists():
        print("Output file exists and will be overwritten")
        # logging.info("Output file exists and will be overwritten")
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
    # logging.debug(data)
    return data



output_fn = file_validation(snakemake.output[0])
df = process_bam(
    bam=snakemake.input[0],
    bed=None,
    stout=False,
    threads=snakemake.threads,
    tol=snakemake.params.tol
)

if isinstance(df, pd.DataFrame):
    # logging.info(f"Writting to {args.output}")
    df.to_csv(output_fn, index=False)
# logging.info("Finished!")
