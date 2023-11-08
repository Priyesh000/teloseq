#!/usr/bin/env python 
        
        
import pysam
import pyfastx
from pyfastx import Fasta 
from collections import defaultdict
import click
from pathlib import Path


class WriteFasta:
    '''Write each arms to sep file '''
    fa = {}
    read_count = defaultdict(int)
    def __init__(self, filename, output_dir, mode='w'):
        self.output_dir = Path(output_dir)
        if not self.output_dir.exists():
            self.output_dir.mkdir(parents=True)
        self.filename = str(self.output_dir / f'{filename}.{{}}.fasta')
        self.mode = mode 
    
    def __enter__(self):
        return 
    
    def __exit__(self, type, value, traceback):
        '''close open files'''
        for chromo, f in self.fa.items():
            print(f'Closing: {f.name}')
            f.close()
            ## if too few reads than delete
            if self.read_count[chromo] < 5:
                Path(f.name).unlink() ## delete file 
        self.clear()

    def write(self, chrom, name, seq):
        if not chrom in self.fa:
            self.fa[chrom] = open(self.filename.format(chrom), mode=self.mode)
            
        # elif self.mode == 'a':
        #     self.fa[chrom] = open(self.filename.format(chrom), mode=self.mode)
        elif self.fa[chrom].closed:
            raise IOError('file is closed')
        self.fa[chrom].write(f'>{name}\n{seq}\n')
        self.read_count[chrom] += 1
    
    def clear(self):
        self.fa = {}
        self.read_count = defaultdict(int)
    

@click.command()
@click.argument('fastafile')
@click.argument('bamfile')
@click.argument('output')
@click.argument('prefix')
@click.option('--skip_supplementary', is_flag=True,  help='skip supplementary alignments [default: True]')
@click.option('--skip_secondary', is_flag=True, help='skip secondary alignments [default: True]')
@click.option('--min_mq', type=int, default=40, help='Filter alignment less than X [default: 40]')
@click.option('--tol', type=int, default=20_000, help='Consider reads N bp from the ends [default: 20_000]')
def main(fastafile, bamfile, output, prefix, skip_supplementary, skip_secondary, min_mq, tol):
    
    print(f'pyfastx version: {pyfastx.version}')
    FA = Fasta(fastafile, build_index=True)
    ## bedtools getfasta rename fasta header to [read_id:corrds ] 
    true_names = {x.split(':')[0]: x for x in FA.keys()}

    # print(FA.keys())
    output = Path(output)

    if output.exists():
        output.unlink()  # delete if exists
        output.mkdir(parents=True)
    else:
        output.mkdir(parents=True)


    samfile = pysam.AlignmentFile(bamfile, 'rb')

    write2fasta = WriteFasta(prefix, output, 'w')
    with write2fasta:
        for rec in samfile:
            if rec.is_unmapped:
                continue
            if rec.is_supplementary and skip_supplementary:
                continue
            if rec.is_secondary and skip_secondary:
                continue
            ## mapping qulity filter
            if rec.mapping_quality <= min_mq:
                continue

            read_id = rec.query_name.split(':')[0]

            refname = rec.reference_name
            reflen = samfile.get_reference_length(refname)
            rstart = rec.reference_start
            rend = rec.reference_end

            if min(rstart, rend) <= tol:
                arms = 'P'
            elif max(rstart, rend) >= (reflen - tol):
                arms = "Q"
            else:
                continue
            refname = f'{refname}_{arms}'
            if not read_id in true_names.keys():
                # print(f'{read_id} not in fasta/telomeric')
                continue
            seq = FA[true_names[read_id]]
            ## write to fasta
            write2fasta.write(refname, read_id, seq)

    ## Close bam
    samfile.close()

if __name__ == "__main__":
    main()
