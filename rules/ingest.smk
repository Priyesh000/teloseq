import numpy as np


rule import_fastq_and_filter:
    input: lookup_value('fastq_path', basecalls_df)
    output: paths.basecalls.fastq 
    params: 
        min_length = config['FILTER']['min_length'],
        max_length = config['FILTER']['max_length']
    conda: 'envs/fastcat.yml'
    shell:
        'fastcat --min_length {params.min_length} --max_length {params.max_length} {input} | gzip -1 > {output}'
    

rule import_tsv:
    input: lookup_value('anchor_tsv', ref_df)
    output: paths.refgenome.tsv
    shell:
        'cp {input} {output}'

# rule tsv_bed:
#     input: paths.refgenome.tsv
#     output: paths.refgenome.bed
#     shell:
#         'awk \'{{print $1"\t"$2"\t"$3}}\' {input} > {output}'



    # conda: 'envs/nanofilt.yml'
    # shell:
    #     '''
    #     if [[ "{input}" == *.gz ]];then  
    #         gzip -dc {input} | NanoFilt --length {params.min_length} --maxlength {params.max_length} | gzip -1 > {output}
    #     else
    #         NanoFilt --length {params.min_length} --maxlength {params.max_length} {input} | gzip -1 > {output}
    #     fi
    #     '''


rule fq2fa:
    input: paths.basecalls.fastq
    output: paths.basecalls.filter
    conda: 'envs/seqkit.yml'
    threads: 4
    shell:
        'seqkit fq2fa -j {threads} -o {output} {input} '

rule virtual_digest:
    input: paths.basecalls.filter
    output: paths.basecalls.fasta
    params:
        enzyme = 'None' if pd.isnull(lookup_value('enzyme', basecalls_df)) else lookup_value('enzyme', basecalls_df),
        out = lambda w, output: Path(output[0]).parent,
        prefix = lambda w, output: Path(output[0]).name.replace('.digest.fasta',''),
    threads: 4
    shell:
        'python scripts/Restriction_digest_reads.py ' 
        '--enzyme {params.enzyme} '
        '--output {params.out} '
        '--prefix {params.prefix} '
        '{input} '



rule import_genomes:
    input: lookup_value('refgenome_path', ref_df)
    output: paths.refgenome.fasta
    shell:
        '''
        if [[ "{input}" == *.gz ]];then  
            gzip -dc {input} > {output} 
        else
            cp {input} {output}
        fi
        '''

# def get_seq_summary():
#     p = lookup_value('summary_path', basecalls_df)        
#     print(str(p))

#     if p is [None, '', np.nan] or not Path(p).exists():
#         return None
#     else:
#         return str(p)

# rule import_seq_summary:
#     output: paths.basecalls.summary
#     params:
#         fn = lookup_value('summary_path', basecalls_df)
#         #fn = get_seq_summary()
#     #conda: 'env/seqkit.yml'
#     run:
#         if params.fn is None:
#             shell('echo "read_id,sequence_length_template,mean_qscore_template"| > {output}.tmp')
#             #shell('seqkit {input} | sed -e "s/\s\+/,/g" | gzip -1 > {output}')
#             shell('seqkit  {input} | sed -e "s/\s\+/,/g" >> {output}.tmp')
#             shell('cat {output}.tmp | gzip -1 > {output}; rm {output}.tmp')
#         else:
#             shell('cp {params.fn} {output}')

rule per_read_stats:
    output: paths.basecalls.stats
    input: paths.basecalls.fastq
    threads: 5
    shell: 
        'seqkit fx2tab -j {threads} --header-line --name --length --only-id --avg-qual {input} | sed -e "s/\s\+/,/g" >> {output}'


rule get_cytoBands:
    output: paths.misc.chromo_arms 
    run:
        fn = Path('cytoBand.txt.gz')
        if not fn.exists():
            shell('wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz')
        bands = pd.read_csv(fn, sep='\t',header=None, usecols=range(4), names=['chromo','start','end', 'cband'])
        bands['arm'] = bands['cband'].str.extract('([p|q])',expand=True)
        bands['chromo_arms'] = bands.apply(lambda x: f"{x['chromo']}.{x['arm']}", axis=1)
        arms_df = (bands 
           .groupby(['chromo', 'chromo_arms'])
           .agg({'start':min, 'end':max})
           .reset_index()[["chromo","start","end", "chromo_arms"]]
           )
        arms_df.to_csv(output[0], index=False, sep='\t', header=None)


# rule add_chromo_arms:
#     input:
#         arms = paths.misc.chromo_arms,
#         stats = paths.mapping.wmap_stat
#     output: paths.mapping.wmap_stat_arms
#     run:
#         from pybedtools import BedTool
#         import pandas as pd 

#         def reorder_list(list1, list2):
#             tmp = [*list1]
#             for item in list2:
#                 if item not in list1:
#                     tmp.append(item)
#             return tmp

#         arms_bed = BedTool(input.arms)
#         bt_col = ['rname', 'rstart', 'rend', 'read_id']
#         df = pd.read_csv(input.stats)
#         new_col = reorder_list(bt_col, df.columns.values)
#         df = df[new_col]  ## set columns order 
#         bam_bt = BedTool.from_dataframe(df)
#         res = arms_bed.intersect(bam_bt, wb=True )
#         res_df = res.to_dataframe(disable_auto_names=True, header=None)
#         #print(type(result), result.columns, dir(result))
#         res_df= res_df.iloc[:,3:]
#         res_df.columns = ['arms'] + new_col
#         res_df.to_csv(output[0], index=False)

# rule LRA_add_chromo_arms:
#     input:
#         arms = paths.misc.chromo_arms,
#         stats = paths.mapping.lra_stat
#     output: paths.mapping.lra_stat_arms
#     run:
#         from pybedtools import BedTool
#         import pandas as pd 

#         def reorder_list(list1, list2):
#             tmp = [*list1]
#             for item in list2:
#                     if item not in list1:
#                         tmp.append(item)
#             return tmp
            
#         arms_bed = BedTool(input.arms)
#         bt_col = ['rname', 'rstart', 'rend', 'read_id']
#         df = pd.read_csv(input.stats)
#         new_col = reorder_list(bt_col, df.columns.values)
#         df = df[new_col]  ## set columns order 
#         bam_bt = BedTool.from_dataframe(df)
#         res = arms_bed.intersect(bam_bt, wb=True )
#         res_df = res.to_dataframe(disable_auto_names=True, header=None)
#         res_df= res_df.iloc[:,3:]
#         res_df.columns = ['arms'] + new_col
#         res_df.to_csv(output[0], index=False)

# rule refgenome_restricted:
#     input: 
#         ref =  paths.refgenome.fasta,
#         bed = paths.refgenome.restricted
#     output: paths.refgenome.restricted_fasta
#     conda: 'envs/bedtools.yml'
#     shell:
#         'bedtools getfasta -fi {input.ref} -bed {input.bed} -fo {output}'

rule agg_basecall_stats:
    input: 
        stats = paths.basecalls.stats,
        telo = paths.telomeric.lor
    output: paths.basecalls.agg_stats 
    run:
        import pandas as pd
        import numpy as np

        ## Fastest
        def N50(a, perc=0.5):
            a = np.array(a)
            a[::-1].sort() ## sort in deascending order
            csum = a.cumsum()
            total = csum.max() 
            nx_idx = np.where(csum == csum[csum>=(total*perc)].min())
            return a[nx_idx][0]
        
        #id,length,avg.qual
        df = pd.read_csv(input.stats)
        df = df.drop_duplicates(subset='#id')
        telomere = pd.read_csv(input.telo, header=None, names=['read_id'])
        telomere['clean_read_id'] = telomere['read_id'].apply(lambda x: x.split('/')[0])
        df.loc[:,'is_telomere'] = False
        df.loc[df['#id'].isin(telomere.clean_read_id), 'is_telomere'] = True
        def agg_stats(x) -> pd.Series:
            d = {
                'total_reads': len(x),
                'total_bases_mb': x['length'].sum()/1e6,
                'read_min_len': x['length'].min(),
                'read_mean_len': x['length'].mean(),
                'read_median_len': x['length'].mean(),
                'read_n50': N50(x['length']),
                'read_max_len': x['length'].max(),
            }
            return pd.Series(d)

        agg_df = df.groupby('is_telomere').apply(agg_stats)
        agg_df.reset_index().to_csv(output[0], index=False)

        
    