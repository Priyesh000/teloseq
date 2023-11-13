import argparse
import re
from collections import Counter, defaultdict
import pandas as pd
from pysam import FastxFile
from typing import Tuple, Dict
from tqdm import tqdm


def generate_variant_pattern(kmer: str) -> str:
    """Generate regex pattern for variant repeats"""
    bases = {
        'N': 'ATGC',
        'B': 'CGT',
        'D': 'AGT',
        'H': 'ACG',
    }
    pattern = [f'[{bases[base]}]' if base in bases.keys() else base for base in kmer]
    # print(pattern)
    return f"({''.join(pattern)})"


def count_kmers(seq: str, ksize: int=6) -> dict:
    """Count kmers in a sequence"""
    kmers = {}
    for i in range(len(seq) - ksize + 1):
        kmer = seq[i:i+ksize]
        if kmer in kmers:
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1
    return kmers


def reduce_to_common_kmer(kmer: str) -> str:
    """Reduce variant repeats to common kmer"""
    return max([kmer[n:]+kmer[:n] for n in range(len(kmer))])
    

def reduce_kmer_counts(kmer_counts: dict) -> dict:
    """Reduce variant repeats to common kmer"""
    reduced_kmer_counts = defaultdict(int)
    for kmer, count in kmer_counts.items():
        reduced_kmer = reduce_to_common_kmer(kmer)
        reduced_kmer_counts[reduced_kmer] += count
    return reduced_kmer_counts


def reverse_complement(kmer:str, comp: bool=True) -> str:
    """Reverse complement a kmer"""
    comp_bases = dict(list(zip('ATCGN', 'TAGCN')))
    new_kmer = []
    for base in kmer[::-1]:
        if comp:
            base = comp_bases[base]
        new_kmer.append(base)
    return ''.join(new_kmer)


def find_motifs(pattern: str, seq: str) -> dict:
    """Find motifs in a sequence"""
    return dict(Counter(re.findall(f'{pattern}', seq, flags=re.I)))


def read_bed_file(bedfile: str) -> Dict[str, Tuple[int, int, str]]:
    '''
    Bed file contrains the read_id, telomere start, telomere end, motif.
    Read bed file format and return a dict of: 
    read_id, start, end, motif[TTAGGG|CCCTAA]
    '''
    store = {}
    with open(bedfile) as fh:
        for record in fh:
            read_id, start, end, motif = record.strip().split('\t')[:4]
            store[read_id] = (int(start), int(end), motif)
    return store


def main():

    coords = read_bed_file(args.bed)
    data = []
    for rec in tqdm(FastxFile(args.fasta)):
        read_id = rec.name
        seq = rec.sequence
        if read_id in coords:
            telomere_start, telomere_end, motif = coords[read_id]
            telomere = seq[telomere_start:telomere_end + 1]

            if motif not in ['CCCTAA', 'TTAGGG']:
                raise ValueError(
                    f'Unknown motif: {motif}, only CCCTAA and TTAGGG are supported')

            if motif == "CCCTAA":
                telomere = reverse_complement(telomere, comp=True)
            
            motifs_counts = count_kmers(telomere, ksize=6)
            motifs_counts = reduce_kmer_counts(motifs_counts)
            # print(motifs_counts)
            # pattern = generate_variant_pattern(args.search_motif)
            # motifs_counts = find_motifs(pattern, telomere)

            rec_dict = {
                'read_id': read_id,
                'start': telomere_start,
                'end': telomere_end,
                'motif': motif,
                'telomeric_len': len(telomere),
                'read_len': len(seq),
            }
            rec_dict.update(motifs_counts)
            data.append(rec_dict)
            # exit()
    df = pd.DataFrame(data)
    df.to_csv(args.output, index=False)


if __name__ == "__main__":

    # assert reduce_to_common_kmer("TAGGGT") == "TTAGGG"
    parser = argparse.ArgumentParser(description='Count telomeric variants')
    parser.add_argument('fasta', metavar='File', help='Fasta reads')
    parser.add_argument('bed', metavar='File',
                        help='Bed file with telomere coords')
    parser.add_argument('output', metavar='File', help='output file')
    parser.add_argument('--search_motif', default='NNNNN',
                        type=str, metavar='String', help='default[NNNGGG]')

    args = parser.parse_args()
    main()
