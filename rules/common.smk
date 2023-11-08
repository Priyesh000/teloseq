from box import Box
import yaml
from pathlib import Path



def create_input_dataframes():
    """Create input dataframes to be igested by snakemake"""
    def read_csv(path, index=None):
        df = pd.read_csv(path, sep=',', comment='#')
        if index:
            df = df.set_index(index, drop=False)
        return df  

    ref_df = read_csv(config['GENOMES'], 'refgenome_id')
    basecalls_df = read_csv(config['BASECALLS'])
    enzymes_df = read_csv(config['ENZYMES']) ## Enzymes list
    basecalls_df = basecalls_df.merge(enzymes_df, on='run_id', how='left')
    basecalls_df = basecalls_df.set_index(['run_id','dataset_id'], drop=False)

    ref_l = ref_df.refgenome_id.unique()
    basecall_ref_l = basecalls_df.refgenome_id.unique()
    missing = []
    for ref in  basecall_ref_l:
        if ref not in basecall_ref_l:
            missing.append(ref)
    if len(missing)>0:
        raise ValueError(f"The following, {missing} refgenome_id from {config['GENOMES']} does not match The refgenome_id in {config['BASECALLS']}")

    ## Expand multi ref
    multi_ref = (
        basecalls_df
        .refgenome_id
        .str
        .split(';', expand=True)
        .stack(dropna=True)
        .rename('ref_ids')
        .to_frame()
        .droplevel(2)

    )
    ## Merge multi ref with basecalls_df
    basecalls_df_stacked = (
        basecalls_df
        .drop('refgenome_id',axis=1)
        .join(multi_ref)
        ).rename(columns={'ref_ids': 'refgenome_id'})
    
    # print(basecalls_df_stacked)
    # print('## ================================== ##')
    # print()
    ## Merge ref_df and basecalls_df
    mapping_df = (
        pd.merge(
            basecalls_df_stacked,
            ref_df.reset_index(drop=True),
            on='refgenome_id',
            how='inner'
        )
    )

    if len(mapping_df) != len(basecalls_df_stacked):
        print(mapping_df)
        raise ValueError(f"The refgenome_id from {config['GENOMES']} does not match The refgenome_id in {config['BASECALLS']}")
    mapping_df = mapping_df.set_index(['run_id','dataset_id', 'refgenome_id'], drop=False)  
    # print(basecalls_df, ref_df, mapping_df)
    return basecalls_df, ref_df, mapping_df

def create_path_accessor(file_layout, prefix: Path) -> Box:
    """Create a Box to provide '.' access to heirarchy of paths"""
    data = yaml.load(Path(file_layout).open(), Loader=yaml.SafeLoader)
    paths = {}
    for directory in data.keys():
        paths[directory] = {}
        for file_alias, file_name in data[directory].items():
            p = str(prefix / directory / file_name)
            if '{{kmer}}' in  p:
                p = p.replace('{{kmer}}', f'Kmer{config["KMC"]["kmer_count"]["k"]}') 
            #if '{refgenome_id}' in p:
            #    p = p.replace('{refgenome_id}', config['refgenome']['refgenome_id'])
            paths[directory][file_alias] = str(p)
    return Box(paths, frozen_box=True)


def to_log(path: str) -> str:
    """Log file location based on output file"""
    return str(outdir / "logs" / path) + ".log"

def to_benchmark(path: str) -> str:
    """Log file location based on output file"""
    return str(outdir / "benchmarks" / path) + ".bench.txt"

def to_prefix(path: str, components=2) -> str:
    """Strip trailing extensions to create an output prefix"""
    return path.rsplit('.', components)[0]


def expand_rows(path: str, df: pd.DataFrame):
    """Expand a templated path string with values from a dataframe"""
    res = df.apply(lambda x: path.format(**x.to_dict()), axis=1)
    return list(res)

def lookup_value(column, df, wildcard_index=None):
    """Use wildcards to 'lookup' a value in a dataframe. The wildcard keys must
    match the dataframe index.
    """
    index_names = tuple(df.index.names)
    if not column in df.columns:
        raise ValueError(f'Column: {column} not in {df.columns}')
    
    def _inner(wildcards):
        if wildcard_index and len(index_names) == 1:
            
            return df.loc[wildcards[wildcard_index], column].values[0]
        if len(index_names) == 1:
            return df.loc[wildcards, column].values[0]
        else:
            ## pd.xs requires a multi index
            row = df.xs(
                tuple(wildcards[k] for k in index_names),
                level=index_names,
                drop_level=True)

            if (len(row) > 1):
                raise ValueError(f'More than one enter {row}')
            return row[column].values[0]
    
    return _inner


def get_opts( config, flag='--',fill_gap=' ', skip_opts=['threads', 't'] ):
    opts = []
    flag_in = flag
    def _mix_flag(flag, key, word_len=1):
        if flag == 'mix':
            if len(key) > word_len:
                return '--'
            else:
                return '-'
        else:
            return flag

    for key, val in config.items():
        if key in skip_opts:
            continue
        if isinstance( val, bool ):
            if val:
                flag = _mix_flag(flag_in, key)
                opts.append(f'{flag}{key}')
        else:
            flag = _mix_flag(flag_in, key)
            opts.append( f'{flag}{key}{fill_gap}{val}' )
    return ' '.join( opts ) 

def kmc_opts(config, flag='-', skip=[]):
    '''
    config: Dict
    return: string with -flag and value join ie -ci80 -t20
    '''
    tmp = []
    for k, v in config.items():
        if k in skip:
            continue
        if isinstance(v, bool):
            if v:
                tmp.append('{}{} '.format(flag, k))
        else:
            tmp.append('{}{}{} '.format(flag, k, v))
    return ' '.join(tmp)

def expand_rows_with_lookup(paths, df, is_aggreated=False):
    """
    Use wildcards to 'lookup' a value in a dataframe. The wildcard keys must
    match the dataframe index. 
    Useful when combining multiple files by dropping a wildcard
    """
    index_names = tuple(df.index.names)
    def __inner(wildcards):

        rows = df.xs(
                wildcards,
                level=list(wildcards.keys()),
                drop_level=True)
        # print(rows)
        results = [
            paths.format(**x._asdict()) 
            for x in rows.reset_index(drop=True).itertuples()]
        if not is_aggreated:
            assert(len(rows) == 1)             
        return set(results)

    return __inner

def get_column_value(column, df, wildcards):
    res = df.loc[(wildcards), column]
    if isinstance(res, str):
        return res
    
class BedFmt:
    index = 0
    def __init__(self,chrom, start, end):
        self.chrom = chrom
        self.start = int(start) if start else None 
        self.end = int(end) if end else None
        BedFmt.index +=1
    
    # def __repr__(self):
    #     return f"{self.__class__}: {' '.join(self.__dict__.values())}" 

    @classmethod
    def reset_index(cls):
        cls.index = 0

    @classmethod
    def from_file(cls, records):
        chrom, start, end, *_ = records.split('\t')
        return cls(chrom, start, end)
    
    @classmethod
    def from_region(cls, region):
        chrom, coords = region.split(':')
        start, end = coords.split('-')
        return cls(chrom,start, end)

    def to_region(self):
        return f"{self.chrom}:{self.start}-{self.end}"
    
    def to_tab(self):
        return '\t'.join(map(str,[self.chrom, self.start, self.end]))
    
    @property
    def shift_one_base(self):
        self.start += 1
        self.end +=1
        
    @property
    def size(self):
        return self.end - self.start

    def to_dict(self):
        return self.__dict__
    
    def is_overlap(self, target):
        t1,t2 = target
        r1,r2 = self.start, self.end 
        if t1<= r1 <= t2 <= r2:
            return True ## 3' overhange
        elif r1<= t1 <= r2 <= t2:
            return True ## 5' overhange
        elif t1 <= r1<= r2 <= t2:
            return True  ## inside
        elif r1 <= t1 <= t2 <= r2:
            return True
        else:
            return False

        
class ReadBedFile:
    def __init__(self, filename):
        self.filename  = filename
        
    def __enter__(self):
        BedFmt.index = 0  ## reset index
        self.file_obj = open(self.filename)
        for records in self.file_obj:
            yield BedFmt.from_file(records)
    
    def __exit__(self, error_type, value, traceback):
        
        self.file_obj.close()
        print('File closed')
        if error_type:
            raise error_type(value, traceback)
        return True