import pandas as pd
import click

def read_tsv(fn):
    df =  pd.read_csv(fn, sep='\t', 
        usecols=["contig", "telomere_start", 'telomere_end',
                 "subtelomere_start", "subtelomere_end", 'chromosome', 'arm']
    )
    df['chrom_arm'] = df['chromosome'] + "." + df['arm']
    ## Get the start and end of the telomere and subtelomere
    df['start'] = df.apply(lambda x: 
                    x['subtelomere_start'] if x['subtelomere_start'] < x['telomere_start'] else x['telomere_start'], axis=1)
    df['end'] = df.apply(lambda x: 
                    x['subtelomere_end'] if x['subtelomere_end'] > x['telomere_end'] else x['telomere_end'], axis=1)
    # df.drop(columns=['chrom', 'arm'], inplace=True)
    # df['length'] = df['end'] - df['start']
    return df[['contig', 'start', 'end', 'chrom_arm']]

@click.command()
@click.argument('archor_tsv', type=click.Path(exists=True, dir_okay=False))
@click.argument('output_filename', type=click.Path(dir_okay=False))
def main(archor_tsv: str, output_filename: str):
    df = read_tsv(archor_tsv)
    df.to_csv(output_filename, sep='\t', index=False, header=False)

if __name__ == "__main__":
    main()
