from pathlib import Path
import click
from Bio import Restriction
from Bio.Restriction import Analysis, RestrictionBatch
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from textwrap import fill
from tqdm import tqdm
from pysam import FastxFile
import logging 
from sys import stderr
import concurrent.futures

## logging 
logging.basicConfig(
            stream=stderr, 
            format='%(asctime)s %(levelname)s: %(message)s',
            datefmt='%d/%m/%Y %H:%M:%S',
            level=logging.INFO)

def get_enzyme_obj(enzyme: str)-> Restriction:
    """Get Restriction Enzyme object"""
    restriction_enzyme = getattr(Restriction, enzyme, None)
    if restriction_enzyme:
        return restriction_enzyme
    raise ValueError(f'Could not find {enzyme}. Check capitalization! ')

def convert_seq_obj(seq: str)-> Seq:
    """Convert sequence string to Seq object"""
    return Seq(seq, IUPACAmbiguousDNA())
        
def digest_sequence(rec, enzyme_obj):
    seq_obj = convert_seq_obj(rec.sequence)
    ori_len = len(rec.sequence)
    results = enzyme_obj.catalyse(seq_obj, linear=True)
    for n, seq in enumerate(results):
        name = f"{rec.name}/{enzyme_obj.__name__}/{n}"
        data = [
            name, ori_len, len(seq),enzyme_obj.__name__,
            enzyme_obj.elucidate(), len(results)
            ]
        yield name, seq, data ## new name, digested seq, meta-data

def fasta_format(name,seq, width=80):
    return f">{name}\n{fill(str(seq), width=width)}\n"

@click.command()
@click.argument('fasta')
@click.option('--enzyme', '-e', default=None, type=str, help='Enzyme')
@click.option('--output', '-o', default=None, type=str, help='Output directory')
@click.option('--prefix','-p',  default=None, type=str, 
    help='Prefix filename  <*.digest.tsv, *.digest.fasta> default[input filename *.digest*]')
def digest(fasta, enzyme , output, prefix):
    """Digest all reads with a restriction Enzyme if none is given the original file is outputed"""
    if output is None:
        output = Path().cwd()

    output = Path(output)
    if not output.exists():
        output.mdkir()
    if prefix is None:
        prefix = Path(fasta).with_suffix('').name
    datalogs = output/f"{prefix}.digest.tsv"
    fasta_out = output/f"{prefix}.digest.fasta"
    COLOMN_NAMES = ('read_id', 'orignal_len', 'digested_len', 'enzyme', 'site', 'num_fragments')

    with FastxFile(fasta) as fa, open(fasta_out,'w') as fa_out, open(datalogs,'w') as dlog:
        enzyme_obj = None
        if enzyme:
            if enzyme.lower() != "none": ## None as string
                enzyme_obj = get_enzyme_obj(enzyme)
                logging.info(f'Using {enzyme_obj.__name__} for restriction digest')
            else:
                enzyme_obj = None
                logging.info(f'No enzyme given')
        else:
            logging.info(f'No enzyme given')
        dlog.write('\t'.join(COLOMN_NAMES)+'\n')  
        for rec in tqdm(fa, 'Progressing Fasta'):

            if enzyme_obj is None:
                fa_out.write(str(rec)+'\n')
            else:
                ## TODO: Multiple processing, if not slow
                for name, seq, data in digest_sequence(rec, enzyme_obj):

                    dlog.write('\t'.join(map(str,data))+"\n")
                    fa_out.write(fasta_format(name, seq))

if __name__ == "__main__":
    digest()
    
