import pandas as pd
import numpy as np
from pathlib import Path
from typing import Literal, List, Dict, Union, Tuple
import pysam
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from loguru import logger
import click
from tqdm import tqdm
import polars as pl
import pyarrow

## TODO:
# 1. ingest bam
# 2. get telomere arm location from anchored bed file
# 3. get telomere arm sequence
# 4 Generate kmer
#   a. Collapse canonical kmers
#   b. Encode kmers 0 normal, 1 variant
#   c. Reverse complement P arms if needed
# 5. Filter reads length
# 6. Generate figures
# 7. Save dataframe to file

__updated__ = "2023-08-31"


def get_target_segment(pairs: List[Tuple[int, int]], target_start: int, target_end: int) -> List:
    """Get the target segment from read.pairs"""
    pos = []
    for pair in pairs:
        if pair[1] >= target_start and pair[1] <= target_end:
            pos.append((pair))
    return pos


def reverse_complement(seq: str) -> str:
    """Reverse complement sequence"""
    seq_dict = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join([seq_dict[base] for base in reversed(seq)])


def targeted_regions(loc_f: str) -> pd.DataFrame:
    logger.info(f"Loading anchored bed file: {loc_f}")
    df = pd.read_csv(loc_f, sep="\t")
    logger.info(f"Anchored bed file shape: {df.shape}")
    return df


def alpha_max(kmer: str) -> str:
    """Reduce variant repeats to common kmer"""
    return max([kmer[k:] + kmer[:k] for k in range(len(kmer))])


def get_kmers(seq: str, kmersize: int = 6) -> List[str]:
    """Split sequence into kmers of size kmersize and reduce to common kmer"""
    return [alpha_max(seq[i : i + kmersize]) for i in range(len(seq) - kmersize + 1)]


def is_canonical(kmer: str) -> int:
    """Determine if a kmer is canonical or variant"""
    return 0 if alpha_max(kmer) == "TTAGGG" else 1


def encode_kmers(kmers: List[str]) -> np.array:
    """Encode kmers as canonical or variant"""
    return np.array([is_canonical(kmer) for kmer in kmers])


def process_bam(bamfile: Path, anchor_bed: Path, min_seqeunce_length: int = 3000) -> pl.DataFrame:
    """ """
    loc_df = targeted_regions(anchor_bed)
    records = []
    logger.info(f"Loading bam file: {bamfile}")
    with pysam.AlignmentFile(bamfile, "rb") as fh:
        for idx, dfx in loc_df.iterrows():
            if not dfx.contig in fh.references:
                continue
            for read in tqdm(
                fh.fetch(dfx.contig, dfx.telomere_start, dfx.telomere_end),
                desc=f"Processing {dfx.contig}: {dfx.telomere_start}-{dfx.telomere_end}",
            ):
                ## get aligned pairs return query/reference position list(query, reference)
                pairs = read.get_aligned_pairs(matches_only=True)
                segment = get_target_segment(pairs, dfx.telomere_start, dfx.telomere_end)
                seq_segment = read.query_sequence[segment[0][0] : segment[-1][0]]
                if dfx.arm == "P":
                    ## Minimap bam stores read in forward strand orientation.
                    # Therefore, canonal kmer would be CCCTAA
                    seq_segment = reverse_complement(seq_segment)
                records.append(
                    {
                        "read_id": read.query_name,
                        "qstart": segment[0][0],
                        "qend": segment[-1][0],
                        "rstart": segment[0][1],
                        "rend": segment[-1][1],
                        "is_reverse": read.is_reverse,
                        "seqlen": len(seq_segment),
                        "telomere_start": dfx.telomere_start,
                        "telomere_end": dfx.telomere_end,
                        "chromosome": dfx.chromosome,
                        "arm": dfx.arm,
                        "seq": seq_segment[:min_seqeunce_length],
                    }
                )
    logger.info("Building dataframe")

    df = pl.DataFrame(records)
    logger.info(f"Dataframe shape: {df.shape}")
    logger.info("Encoding kmers")
    # df = df.with_columns(
    #     (pl.col('seq')
    #         .apply(lambda x: encode_kmers(get_kmers(x, 6)))
    #         .alias('encode')
    #     )
    #     )
    ## .apply is deprecated. Use .map_elements instead
    df = df.with_columns(
        (pl.col("seq").map_elements(lambda x: encode_kmers(get_kmers(x, 6))).alias("encode"))
    )
    return df
    # return pd.DataFrame(records)


# def tvr_sig(dfx: pd.DataFrame) -> float:
#     """Normalize TVR signal"""
#     logger.info(f'Calculating TVR signal for {dfx.chromosome} {dfx.arm}')
#     return round(100 * (dfx.sum()/dfx.shape[0]),3) ## normalize by number of reads


def tvr_sig(encoding) -> pl.DataFrame:
    """Normalize TVR signal"""
    total_sum = np.stack(encoding).sum(0)
    normalized_sum = 100 * (total_sum / len(encoding))
    return normalized_sum


def plot_tvr_signal(x: str, **kwargs) -> plt.gca():
    """Plot TVR signal"""
    # logger.info(f'Plotting TVR signal {kwargs.get("label")}')
    ax = plt.gca()
    data = kwargs.pop("data", None)
    ax.plot(data[x].values[0])
    return ax


def plot_facet_grid(data: pl.DataFrame, output: str) -> sns.FacetGrid:
    """Plot facet grid. Did not make it generic intentionally to avoid confusion"""
    logger.info("Plotting facet grid")
    chrom_order = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]
    arm_order = ["P", "Q"]
    g = sns.FacetGrid(
        data,
        col="arm",
        row="chromosome",
        aspect=2,
        row_order=chrom_order,
        col_order=arm_order,
    )
    g.map_dataframe(plot_tvr_signal, "signal")
    g.figure.savefig(output, dpi=300)
    return g


def normalize_encoding(
    df: pl.DataFrame, save_raw_data: Path, min_seqeunce_length: int = 3000, min_reads: int = 10
):
    """Normalize TVR signal"""
    ## filter for min sequence length
    df = df.filter(pl.col("seqlen") > min_seqeunce_length)
    ## Read counts after legnth filter
    df = df.with_columns(pl.col("read_id").count().over(["chromosome", "arm"]).alias("read_counts"))
    # df.write_parquet(save_raw_data)
    # logger.info(f"Raw writen to file after length filtering but no read count filter: {save_raw_data}")
    ## filter for min number of reads
    df = df.filter(pl.col("read_counts") > min_reads)

    ## group by chromosome and arm
    norm = df.groupby(["chromosome", "arm"]).agg(
        [
            pl.col("encode").apply(tvr_sig).alias("signal"),
            pl.col("read_id").count().alias("read_counts"),
        ]
    )
    norm = norm.to_numpy()  ## convert to numpy array cos to_pandas is not working
    norm = pd.DataFrame(norm, columns=["chromosome", "arm", "signal", "read_counts"])
    return norm


@click.command()
@click.option(
    "--bamfile",
    "-b",
    type=click.Path(file_okay=True, exists=True, readable=True, dir_okay=False),
    required=True,
    help="Bam file",
)
@click.option(
    "--anchor_bed",
    "-a",
    type=click.Path(file_okay=True, exists=True, readable=True, dir_okay=False),
    required=True,
    help="Anchor bed file",
)
@click.option(
    "--output_dir",
    "-o",
    type=click.Path(file_okay=False, dir_okay=True),
    required=True,
    help="Output directory path",
)
@click.option(
    "--prefix",
    "-p",
    type=str,
    required=True,
    help="Prefix for output file [*.tvr_signal.{png,csv}]",
)
@click.option(
    "--log", "-l", type=click.Path(file_okay=True, dir_okay=False), required=False, help="Log file"
)
@click.option(
    "--min_seqeunce_length",
    "-m",
    type=int,
    default=3000,
    help="Minimum sequence length to consider [default: 3000]",
)
@click.option(
    "--min_reads",
    "-r",
    type=int,
    default=10,
    help="Minimum number of reads per chromosome arm [default: 10]",
)
@click.option(
    "--fmt",
    "-f",
    type=click.Choice(["png", "csv"]),
    default="png",
    help="Output format for figure: \{svg, png\} [default: png]",
)
def main(bamfile, anchor_bed, output_dir, prefix, log, min_seqeunce_length, min_reads, fmt):
    if log:
        logger.add(log, level="INFO")
    output = Path(output_dir, prefix + ".tvr_signal")

    raw_data = str(output) + ".tvr_signal.raw_data.parquet"

    df = process_bam(bamfile, anchor_bed)
    data = normalize_encoding(
        df, save_raw_data=raw_data, min_seqeunce_length=min_seqeunce_length, min_reads=min_reads
    )
    logger.info(f"shape of dataframe: {df.shape}")

    # data = data.to_numpy()
    # data = pd.DataFrame(data, columns=['chromosome', 'arm', 'signal'])
    logger.info(f"Dataframe shape: {data.shape}")
    plot_facet_grid(data, str(output) + f".norm.{fmt}")
    logger.info(f'save figure to {str(output) + f".norm.{fmt}"}')
    data.to_parquet(str(output) + ".norm.parquet", engine="pyarrow")
    logger.info(f'save data to {str(output) + ".norm.parquet"}')


if __name__ == "__main__":
    main()
