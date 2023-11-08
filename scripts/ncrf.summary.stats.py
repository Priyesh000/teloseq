
import pandas as pd
import sys

def rc(seq):
    bases = dict(zip(list('ATGC'), list('TACG')))
    rev  = []
    for base in seq[::-1]:
        rev.append(bases[base])
    return ''.join(rev)

def agg_func(x):

    names = {
        'mstart': x['start'].min(),
        'mend': x['end'].max(),
        'telomere_len': (x['end'].max() - x['start'].min()),
        'motif': rc(x['motif'].max()) if x['strand'].max() == '-' else x['motif'].max(),
        'subtelomere_len': (x['seqLen'].max() - (x['end'].max() - x['start'].min()) ),
        'read_len': x["seqLen"].max(),
        'segments': len(x)

    }
    for i in ['m','mm','i', 'd']:
        names[i] = x[i].sum()

    return pd.Series(names, index=names.keys())

def agg_mean(x):
    names = {
        'read_count': len(x),
        'min_telomere_len' : x['telomere_len'].min(),
        'mean_telomere_len': x['telomere_len'].mean(),
        'max_telomere_len': x['telomere_len'].max(),

        'min_subtelomere_len' : x['subtelomere_len'].min(),
        'mean_subtelomere_len': x['subtelomere_len'].mean(),
        'max_subtelomere_len': x['subtelomere_len'].max(),

        'min_read_len': x['read_len'].min(),
        'mean_read_len': x['read_len'].mean(),
        'max_read_len': x['read_len'].max(),

        'mean_ident': x['ident'].mean()
        }

    return pd.Series(names, index=names.keys())

def read_start_correctly(row: pd.DataFrame, padding: int=200) -> bool:
    """Check if the read starts with correct the motif """
    ## 5' CCCTAA at start of read: 5' <ONT_seq>CCCTAA====> [CENT] 5'<ONT_seq>====TTAGGG> 3'> 
    ## 5 TTAGGG at the end of read 3' <TTAGGG====<ONT_seq>5' [CENT] <====CCCTAA<ONT_seq>5'> 
    if row.motif == 'CCCTAA':
        if row.mstart < padding:
            return True
        else:
            return False
    elif row.motif == "TTAGGG":
        if row.mend > (row.read_len-padding):
            return True
        else:
            return False


def telomere_terminal(row: pd.DataFrame, padding: int = 200) -> bool:
    """check if the telomoere mofit start at terminals ends of read"""
    if row.mstart < padding:
        return True
    elif row.mend > (row.read_len - padding):
        return True
    else:
        return False

def main():
    
    df = pd.read_csv(
        sys.argv[1],
        sep='\t',
        usecols=['seq', 'start', 'end', 'seqLen','motif', 'strand', 'm', 'mm', 'i', 'd']
        )
    
    df = df.rename(columns={'seq': 'read_id'})
    df['motif'] = df.motif.str.replace('telomeric', 'TTAGGG')

    df1 = df.groupby('read_id').apply(agg_func)
    df1['ident'] = 100*(df1.m / df1[['m', 'mm', 'i', 'd']].sum(1))
    
    ## read classification 
    ## >60bp because the shortest digested sub-telomeric len is 63
    df1['has_subtelomere'] = df1.subtelomere_len.apply(lambda x: True if x>500 else False)
    ## add paddding of 200bp for sequence adaptor  
    df1['started_correctly'] = df1.apply(read_start_correctly, axis=1)
    ## Check if alignment starts or ends at the chromomsome terminal ends
    df1['is_terminal'] = df1.apply(telomere_terminal, axis=1)
    
    df1.reset_index().to_csv(sys.argv[2], index=False)

    ## filter: 
    ## segment < 1
    ##  telomere_len > 
    df2 = df1.query('telomere_len >= 1000').groupby('motif').apply(agg_mean)
    print(df2)

if __name__ == "__main__":
    main()
