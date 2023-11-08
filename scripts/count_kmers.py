# from os import W_OK
# import sys
# import subprocess
from pathlib import Path
from typing import Pattern, Union
from numpy.lib.arraysetops import unique
from pysam import FastxFile

import numpy as np
import pandas as pd
from subprocess import check_output
import click
import re

from functools import lru_cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
from functools import partial

from concurrent.futures import ThreadPoolExecutor, as_completed

import regex

REPORT_COLUMNS = [
    "monomer", "motif", "length", "score", "fraction_explained",
    "p", "p_adjusted",
]
REPORT_COLUMNS_ESCAPED = ["#"+REPORT_COLUMNS[0]] + REPORT_COLUMNS[1:]

## progress bar format 
progressbar_template = partial(
    tqdm, bar_format=(
        "{desc}{percentage:3.0f}% ({n_fmt}/{total_fmt}), " +
        "{elapsed}<{remaining}, {rate_fmt}"
    )
)

@click.group()


def cli():
    pass


def run_cmd(cmd):
    print('Executing: ', cmd)
    return check_output(cmd.split())


@cli.command(name='run_jellyfish')
@click.argument('fasta', type=click.Path(exists=True))
@click.argument('prefix', type=str)
@click.argument('wrkdir', type=click.Path(dir_okay=True))
@click.option('--min_k', default=4, type=int, help='Minimum kmer size [default: 4] ')
@click.option('--max_k', default=16, type=int, help='Maxmum kmer size [default: 16]')
@click.option('--threads', default=4, type=int, help='threads [default: 4]')
def jellyfish_count_kmer(fasta: str, prefix: str,  wrkdir: str, min_k: int = 2, max_k: int = 12,
                            threads: int=4, hashsize: str='2G', min_repeat_size: int=2 ):
    """Jellyfish to count kmers. Will create a master.tsv.gz that can be used in \'report\' step"""

    wrkdir = Path(wrkdir)
    if not wrkdir.exists():
        wrkdir.mkdir(parents=True)
    data = []

    for ksize in range(min_k, max_k):
        jf_out = str(wrkdir / f'{prefix}.{ksize}.jf')
        report = str(wrkdir / f'{prefix}.{ksize}.tsv')

        doublets = ksize * min_repeat_size  # Looking for repeating kmer e.g. TTAGGGTTAGGG
        ##Count kmer
        cmd_jf_count = f'jellyfish count -t {threads} -s {hashsize}  -L 0 -C -m {doublets} -o {jf_out} {str(fasta)}'

        ## export to plan text
        cmd_jf_dump = f'jellyfish dump -c -t -L 0 -o {report} {jf_out}'

        cmd_remove_db = f'rm -fr {jf_out}'

        for cmd in [cmd_jf_count, cmd_jf_dump, cmd_remove_db]:
            run_cmd(cmd)
        
        df = pd.read_csv(report, sep='\t', header = None, names=['kmer', 'count'])
        if df.size == 0:
            continue 
        else:
            mask = df.kmer.apply(
                lambda x: True if x[:ksize] * min_repeat_size == x else False)
            df = df[mask]
            if df.size == 0:
                continue
            df['kmer'] = df.kmer.apply(lambda x: x[:ksize])
            # df.loc[:,'kmersize'] = int(ksize)
            df.loc[:, 'length'] = ksize
            data.append(df)
        master_df =pd.concat(data, ignore_index=True)
        
        master_df.to_csv(wrkdir/f'{prefix}.jf.count_filter.tsv.gz', index=False, sep='\t')
    return master_df

## sourced: https://github.com/LankyCyril/edgecas
def validate_motif(motif):
    """Make sure motif conforms to what we've implemented"""
    return (re.search(r'^[ACGT\[\]\.]*$', motif, flags=re.IGNORECASE) is not None)


@lru_cache(maxsize=None)
def revcomp(sequence, ignorecase=True):
    """Reverse-complement a sequence"""
    bases = 'AGCT'
    COMPLEMENTS = dict(zip(bases, bases[::-1]))
    COMPLEMENT_PATTERN = re.compile(r'|'.join(COMPLEMENTS.keys()))
    if ignorecase:
        def matcher(match): return COMPLEMENTS[match.group().upper()]
    else:
        def matcher(match): return COMPLEMENTS[match.group()]
    return COMPLEMENT_PATTERN.sub(matcher, sequence[::-1])



@lru_cache(maxsize=None)
def lowest_common_kmer(kmer):
    """Get ordered lowest common of kmer (e.g., for TTAGGG ==> AGGGTT)"""
    return max(kmer[i:]+kmer[:i] for i in range(len(kmer)))


@lru_cache(maxsize=None)
def common_collapsed_revcomp_kmer(kmer):
    """Get ordered lowest common revcomp or forward of kmer regardless of strand (e.g., for TTAGGG will return AACCCT)"""
    return max(
        lowest_common_kmer(kmer), lowest_common_kmer(revcomp(kmer)),
    )


def custom_order_motif(motif):
    """Custom order motif looks closest to canonical (e.g., for AGGGTTC will return TTCAGGG)
    * AGGGTTC => TTCAGGG
    * TTCAGGG => TTCAGGG
    * CCCTTA => CCCTTA
    * TTACCC => CCCTTA
    """
    a, c, g, t = [motif.count(letter) for letter in "ACGT"]
    is_g = (g > c) or ((g == c) and (t > a))
    if is_g:  ## is G rich
        ## AGGGTTC => TTCAGGG
        if t > 0:
            i = motif.find("T") ## index of first T
            # print('T',i)
        else:
            i = motif.find("G") ## index of first G
            # print('G',i)
    else:
        ## TTACCC => CCCTTA
        if c > 0: ## C rich
            i = motif.find("C") ## index of first C
            # print('C',i)
        else:
            i = motif.find("A") ## index of first A
            # print('A',i)
    if i == -1:
        return min(motif[i:]+motif[:i] for i in range(len(motif)))
    else:
        return motif[i:] + motif[:i]

## https://stackoverflow.com/questions/34947578/how-to-vectorize-fishers-exact-test
## https://towardsdatascience.com/fishers-exact-test-from-scratch-with-python-2b907f29e593
def safe_fisher_exact(count, total_candidate_count, med, total_background_count):
    """Same as running fisher_exact, but returns p=1 if encounters NaNs (i.e., when not enough data to process test)"""
    try:
        return fisher_exact([
            [count, total_candidate_count - count],
            [med, total_background_count - med],
        ])[1]
    except ValueError:
        return 1


def get_motifs_fisher(single_length_report, collapse_reverse_complement=False):
    """Analyze repeat enrichment given the same motif length"""
    lengths = np.unique(single_length_report["length"].values)
    if len(lengths) != 1:
        raise ValueError("`get_motifs_fisher`: multiple lengths found")
    fishery = single_length_report.copy()
    fishery["motif"] = fishery["kmer"].apply(
        common_collapsed_revcomp_kmer if collapse_reverse_complement
        else lowest_common_kmer
    )
    fishery_groupby = fishery[["motif", "count"]].groupby(
        "motif", as_index=False,
    )
    fishery = fishery_groupby.sum()
    iqr = fishery["count"].quantile(.75) - fishery["count"].quantile(.25)
    whisker = fishery["count"].quantile(.75) + 1.5 * iqr
    candidates, background = (
        fishery[fishery["count"] >= whisker].copy(),
        fishery[fishery["count"] < whisker].copy(),
    )
    total_candidate_count, total_background_count = (
        candidates["count"].sum(), background["count"].sum(),
    )
    med = fishery["count"].median()
    candidates["p"] = candidates["count"].apply(
        lambda count: safe_fisher_exact(
            count, total_candidate_count,
            med, total_background_count,
        )
    )
    candidates["length"] = lengths[0]
    return candidates


def analyze_repeats(full_report, collapse_reverse_complement=False, adj="bonferroni"):
    """Analyze repeat enrichment for multiple lengths and apply multiple testing adjustment"""
    candidates = pd.concat([
        get_motifs_fisher(
            full_report[full_report["length"] == length],
            collapse_reverse_complement=collapse_reverse_complement,
        )
        for length in progressbar_template(
            np.unique(full_report["length"].values), unit="k",
            desc="Calculating enrichment ",
        )
    ])
    candidates["p_adjusted"] = multipletests(candidates["p"], method=adj)[1]
    return candidates[
        ["motif", "length", "count", "p", "p_adjusted"]
    ]


def coerce_and_filter_report(analysis, max_p_adjusted):
    """Collapse functionally synonymous entries like TATATA and TATA"""
    motif_mapper, mapped_motifs = {}, set()
    for motif in analysis.sort_values(by="length")["motif"]:
        length_indexer = analysis["length"] > len(motif)
        for larger_motif in analysis.loc[length_indexer, "motif"]:
            if (larger_motif not in mapped_motifs) and (motif in larger_motif):
                if motif*len(larger_motif) == larger_motif*len(motif):
                    if motif not in motif_mapper:
                        motif_mapper[motif] = {motif}
                    motif_mapper[motif].add(larger_motif)
                    mapped_motifs.add(larger_motif)
    synonyms_to_keep = set()
    for synonyms in motif_mapper.values():
        synonym_data = analysis[
            analysis["motif"].isin(synonyms) &
            (analysis["p_adjusted"] < max_p_adjusted)
        ]
        if len(synonym_data):
            synonyms_to_keep.add(
                synonym_data.sort_values(
                    by="count", ascending=False
                ).iloc[0, 0]
            )
    synonyms_to_remove = (
        set.union(set(), *motif_mapper.values()) - synonyms_to_keep
    )
    return analysis[
        (~analysis["motif"].isin(synonyms_to_remove)) &
        (analysis["p_adjusted"] < max_p_adjusted)
    ].copy()


def get_circular_pattern(motif, repeats=2):
    """Convert motif into circular regex pattern (e.g., r'TCGA|CGAT|GATC|ATCG' for TCGA)"""

    atom_pattern = regex.compile(
        r'[ACGT.]|\[[ACGT]+\]', flags=regex.IGNORECASE)
    atoms = atom_pattern.findall(motif)
    if "".join(atoms) != motif:
        raise ValueError("Could not parse motif: {}".format(motif))
    repeated_inversions = {
        "".join(atoms[i:] + atoms[:i]) * repeats
        for i in range(len(atoms))
    }
    pattern = r'|'.join(repeated_inversions)
    RX = regex.compile(pattern, flags=regex.IGNORECASE)
    return RX


def explain_report(filtered_analysis, sequencefile, min_repeats, jobs=1):
    """Calculate fraction of reads explainable by each motif"""
    explained_analysis = filtered_analysis.copy()
    explained_analysis["bases_explained"], total_bases = 0.0, 0
    with FastxFile(sequencefile) as fastx:
        def get_number_of_masked_positions(entry, motifs):
            """get sequence """
            sequence = entry.sequence             
            n_masked_positions_per_motif = {}
            for motif in motifs:
                positions_to_mask = set()

                motifs_pattern = get_circular_pattern(motif, repeats=min_repeats)
                matcher = motifs_pattern.finditer(
                    str(sequence), overlapped=True)
                for match in matcher:
                    positions_to_mask |= set(range(match.start(), match.end()))
                n_masked_positions_per_motif[motif] = len(positions_to_mask)
            return n_masked_positions_per_motif, len(sequence)
        
        with ThreadPoolExecutor(max_workers=1) as pool:
            workers = [
                pool.submit(
                    get_number_of_masked_positions, entry,
                    set(filtered_analysis["motif"]),
                )
                for entry in fastx if entry.sequence
            ]
            iterator = progressbar_template(
                as_completed(workers), total=len(workers),
                desc="Calculating fractions", unit="read",
            )
        
            for worker in iterator:
                n_masked_positions_per_motif, total_seq_bases = worker.result()
                for motif, n_pos in n_masked_positions_per_motif.items():
                    indexer = (
                        explained_analysis["motif"] == motif, "bases_explained",
                    )
                    explained_analysis.loc[indexer] += n_pos
                total_bases += total_seq_bases
    
    return explained_analysis, total_bases


def coerce_to_monomer(motif, min_k):
    """Coerce motif to monomer, e.g. TATA -> TA, CAT -> CAT; this can be used to find functionally synonymous entries too"""
    n = len(motif)
    for i in range(min_k, int(n/2)+1):
        q, r = divmod(n, i)
        if r == 0:
            if motif[:i]*q == motif:
                return motif[:i]
    else:
        return motif


def format_analysis(explained_analysis, min_k, max_motifs, total_bases):
    """Make dataframe prettier"""
    explained_analysis["motif"] = explained_analysis["motif"].apply(
        custom_order_motif,
    )
    explained_analysis["monomer"] = explained_analysis["motif"].apply(
        lambda motif: coerce_to_monomer(motif, min_k=min_k),
    )
    formatted_analysis = explained_analysis.sort_values(
        by=["count", "p_adjusted"], ascending=[False, True],
    )
    formatted_analysis['total_bases'] = total_bases
    formatted_analysis["score"], formatted_analysis["fraction_explained"] = (
        formatted_analysis["count"] / total_bases,
        formatted_analysis["bases_explained"] / total_bases
    )
    formatted_analysis_all = formatted_analysis.copy()
    formatted_analysis = formatted_analysis[REPORT_COLUMNS]   
    formatted_analysis.columns = REPORT_COLUMNS_ESCAPED
    if max_motifs is None:
        return formatted_analysis, formatted_analysis_all
    else:
        return formatted_analysis[:max_motifs], formatted_analysis_all


@cli.command('report')
@click.argument('fasta', type=click.Path(exists=True))
@click.argument('report', type=click.Path(exists=True))
@click.argument('prefix', type=str)
@click.argument('wrkdir', type=click.Path(dir_okay=True))
@click.option('--threads', default=4, type=int, help='threads [default: 4]')
@click.option('--min_k', default=4, type=int, help='Minimum kmer size [default: 4] ')
@click.option('--max_motifs', default=None, type=int,  help='Maximum  number of motifs to report [default: None]')
@click.option('--max_p_adjusted', default=0.05, type=float, help='cutoff adjusted p-value [default: .05]')
@click.option('--min_repeats', default=2, type=int,  help='minimum number of consecutive repeats [default: 2]')
@click.option('--collapse_reverse_complement', is_flag=True, help='collapse reverse complement')
def filter_report(fasta: Path, report: Union[pd.DataFrame, Path], prefix: str, wrkdir: Path, threads: int=4, collapse_reverse_complement: bool = True,
        min_repeats: int = 2, max_p_adjusted: float = 0.05, max_motifs: int=None, min_k: int=4):
    """Filter master jellyfish kmers counts and calculate P value with Fisher exact test"""
    if not isinstance(report, pd.DataFrame):
        if not Path(report).exists():
            raise ValueError(f'{report} do not exists. Please check the path')
        else:
            full_report = pd.read_csv(report, sep='\t')
    else:
        full_report = report


    wrkdir = Path(wrkdir)
    interm = wrkdir / 'interm'
    if not interm.exists():
        interm.mkdir(parents=True)

    outputfile = wrkdir / f'{prefix}.final_report.tsv'
    analysis = analyze_repeats(full_report, collapse_reverse_complement)

    filtered_analysis = coerce_and_filter_report(analysis, max_p_adjusted)

    filtered_analysis.to_csv(interm / f'{prefix}.filtered_analysis.csv.gz')

    explained_analysis, total_bases = explain_report(
        filtered_analysis, fasta, min_repeats, jobs=threads)
    
    explained_analysis.to_csv(interm / f'{prefix}.explained_analysis.csv.gz')

    formatted_analysis, full = format_analysis(
        explained_analysis, min_k, max_motifs, total_bases)
    
    full.to_csv(interm / f'{prefix}.full.csv.gz')

    formatted_analysis.to_csv(outputfile, sep="\t", index=False)
    print(f'written to {outputfile}')



@cli.command('run_full_pipeline')
@click.argument('fasta', type=click.Path(exists=True))
@click.argument('prefix', type=str )
@click.argument('wrkdir', type=click.Path(dir_okay=True))
@click.option('--min_k', default=4, type=int, help='Minimum kmer size [default: 4] ')
@click.option('--max_k', default=16, type=int, help='Maxmum kmer size [default: 16]')
@click.option('--threads', default=4, type=int, help='threads [default: 4]')
@click.option('--max_motifs', default=None, type=int,  help='Maximum  number of motifs to report [default: None]')
@click.option('--max_p_adjusted', default=0.05, type=float, help='cutoff adjusted p-value [default: .05]')
@click.option('--min_repeats', default=2, type=int,  help='minimum number of consecutive repeats [default: 2]')
@click.option('--collapse_reverse_complement', is_flag=True, help='collapse reverse complement')
def count_kmers(fasta, prefix,  wrkdir, min_k, max_k,
                threads,  min_repeats,  max_p_adjusted, max_motifs, collapse_reverse_complement, hashsize='2G'):
    """Run the full pipeline without interaction from use"""
    collapse_reverse_complement = True 
    # max_p_adjusted = 0.05
    # min_repeats = 2
    # max_motifs = None
    wrkdir = Path(wrkdir)
    interm = wrkdir / 'interm'
    if not interm.exists():
        interm.mkdir(parents=True)

    outputfile = wrkdir / f'{prefix}.final_report.tsv'

    full_report = jellyfish_count_kmer(fasta, prefix, wrkdir, min_k, max_k, 
                                       threads, hashsize, min_repeats)
    if full_report is None:
        print(*REPORT_COLUMNS_ESCAPED, sep="\t", file=str(outputfile))

    else:
        filter_report(fasta=fasta,
                      report=full_report, 
                      prefix=prefix, 
                      wrkdir=wrkdir, 
                      threads=threads, 
                      collapse_reverse_complement=collapse_reverse_complement,
                      min_repeats=min_repeats, max_p_adjusted=max_p_adjusted, 
                      max_motifs=max_motifs, min_k=min_k)
        # analysis = analyze_repeats(full_report, collapse_reverse_complement)
        
        # filtered_analysis = coerce_and_filter_report(analysis, max_p_adjusted)

        # filtered_analysis.to_csv(interm / f'{prefix}.filtered_analysis.csv.gz')

        # explained_analysis, total_bases = explain_report(
        #     filtered_analysis, fasta, min_repeats, jobs=threads)

        # formatted_analysis, _ = format_analysis(
        #     explained_analysis, min_k, max_motifs, total_bases)

        # formatted_analysis.to_csv(outputfile, sep="\t", index=False)
        # print(f'written to {outputfile}')

    
if __name__ == "__main__":
    # count_kmers()
    cli()
