import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pybedtools import BedTool
from pathlib import Path    
import click


def munch_methyl_bed(df, min_cov=10):

    colnames = ['cov', 'methylation', 'mod','canconical', 'other_mod', 'delete', 'fail', 'diff', 'nocall']
    df = df.join(df.blockCount.apply(lambda x: pd.to_numeric(pd.Series(x.split(), index=colnames))))
    df.drop('blockCount', axis=1, inplace=True)
    df['chrom_arm'] = df.chrom.apply(lambda x: x.split('/')[0])
    df['arm'] = df.chrom_arm.apply(lambda x: x.split('.')[1])
    df['chromosome'] = df.chrom_arm.apply(lambda x: x.split('.')[0])
    df = df.query('cov >= @min_cov')    
    return df


def plot_methylation(**kwargs):
    data = kwargs.pop("data")
    data = data.sort_values('start')
    ax = plt.gca()
    arm = data.arm.unique()[0]
    if arm == 'Q':
        pos_max = data.end.max() - 15_000
        data = data[data.start > pos_max]
    else:
        data = data[data.start < 15_000]
    ax.plot('start', 'methylation', label='methylation', data=data)
    ax.scatter('start', 'methylation', data=data, s=2)


def get_chrom_ends(chrom_ends_bed, methyl_bed) -> pd.DataFrame:
        
    bed = BedTool(str(chrom_ends_bed))
    meth = BedTool(str(methyl_bed))
    overlaps  = meth.intersect(bed, wa=True)
    return overlaps.to_dataframe() 


def plot_figure(df, savefile_name):
    g = sns.FacetGrid(
        df, 
        row='chromosome', 
        row_order=[f'chr{x}' for x in range(1, 23)],
        col='arm',
        col_order=['P', 'Q'],
        aspect=3,
        sharex=False,
    ).map_dataframe(plot_methylation)

    g.savefig(savefile_name)


@click.command()
@click.argument('bedfile', type=click.Path(file_okay=True, dir_okay=False, exists=True))
@click.argument('methyl_bed', type=click.Path(file_okay=True, dir_okay=False, exists=True))
@click.option('--min-read-cov', '-m', default=10, show_default=True, help='Minimum read coverage')
@click.argument('output_file',)
def main(bedfile, methyl_bed, output_file, min_read_cov):
    """
    
    BEDFILE: 'Path to chrom ends bed file'
    METHTYL_BED: 'Path to methyl bed file from modkit'
    OUTPUT_FILE: 'Save file name'
    """
    df = get_chrom_ends(bedfile, methyl_bed)
    df = munch_methyl_bed(df, min_read_cov)
    df.to_csv(Path(output_file).with_suffix('.csv'), index=False)
    plot_figure(df, output_file)

    

if __name__ == '__main__':
    main()