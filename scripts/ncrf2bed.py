
import argparse
import sys
import pandas as pd 


def reads_ncrf(fn):
    df = pd.read_csv(fn, sep='\t', usecols=['seq', 'start', 'end', 'strand'])
    df = df.rename(columns={'seq': 'read_id', 'strand': 'motif'})
    df['motif'] = df.motif.replace({'+':'TTAGGG', '-':'CCCTAA'})
    return df

def agg_telomere_stats(df):
    def agg_func(x):
        names = {
            "start": x['start'].min(),
            "end":  x['end'].max(),
            'motif': x['motif'].max()
        }
        return pd.Series(names, index=names.keys())
    
    return df.groupby('read_id').apply(agg_func).reset_index()

def main():

    agg = False

    df = reads_ncrf(sys.argv[1])
    if agg:
        agg_telomere_stats(df).to_csv(sys.argv[2],index=False, sep='\t', header=False)
    else:
        df[['read_id', 'start', 'end', 'motif']].to_csv(sys.argv[2],index=False, sep='\t', header=False)

if __name__ == "__main__":
    main()