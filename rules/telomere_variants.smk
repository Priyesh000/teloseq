# from box import Box

### Plan ###
## full length file bam  ==> paths.variants.bam
## combine bams by genome
## extract telomoere from reads
## group reads by chromomsome arms
## count kmers 
## colapse kmer 
## fisher exact test
## plot heat map


rule combine_telomere_only:
    '''Combine telomere sequences into a single file '''
    input: expand_rows_with_lookup(paths.telomeric.telo_only, mapping_df, is_aggreated=True)
    output: paths.telomeric.fa_merged
    shell:
        'cat {input} > {output}'


checkpoint split_by_arms:
    input: 
        fa = paths.telomeric.fa_merged,
        bam = paths.variants.bam
    # output: paths.telomere_variants.fa_arms
    output: directory(str(Path(paths.telomere_variants.fa_arms).parent))
    params:
        prefix = lambda wildcards: Path(paths.telomere_variants.fa_arms).name.replace('.{arms}.fasta','').format(**wildcards) ,
        conf = get_opts(config['KMER_COUNTS']['SPLIT_BY_ARMS'])
    shell:
        'python  scripts/split_telomere_by_arms.py {input.fa} {input.bam} {output} {params.prefix} {params.conf}'


rule run_jellyfish:
    input: 
        fa = paths.telomere_variants.fa_arms 
    output: paths.telomere_variants.jellyfish
    params:
        outdir = lambda x, output: str(Path(output[0]).parent),
        prefix = lambda x, output: str(Path(output[0]).name).replace('.jf.count_filter.tsv.gz', ''),
        conf = get_opts(config['KMER_COUNTS']['TELOMERIC_KMERS'], '--', 
                skip_opts=['t', 'threads', 'max_p_adjusted', 'min_repeats', 'max_motifs'])
    threads: config['KMER_COUNTS']['TELOMERIC_KMERS']['threads']
    shell:
        'python scripts/count_kmers.py run_jellyfish {input.fa} {params.prefix} {params.outdir} {params.conf} --threads {threads}'


rule filter_telomeric_kmers:
    input: 
        fa = paths.telomere_variants.fa_arms,
        counts = paths.telomere_variants.jellyfish 
    output: paths.telomere_variants.report
    params:
        outdir = lambda x, output: str(Path(output[0]).parent),
        prefix = lambda x, output: str(Path(output[0]).name).replace('.final_report.tsv', ''),
        conf = get_opts(config['KMER_COUNTS']['TELOMERIC_KMERS'], '--', 
            skip_opts=['t', 'threads', 'max_k',])
    threads: config['KMER_COUNTS']['TELOMERIC_KMERS']['threads']
    shell:
        'python scripts/count_kmers.py report {input.fa} {input.counts} {params.prefix} {params.outdir} {params.conf} --threads {threads} --collapse_reverse_complement'



def agg_telomeric_kmers_report(wildcards):
    checkpoint_output = checkpoints.split_by_arms.get(**wildcards).output[0] ## Get output dir and fill wildcards
    checkpoint_prefix = checkpoints.split_by_arms.get(**wildcards).rule.params ## Access prefix from rule
    prefix = checkpoint_prefix.prefix(wildcards) ## get prefix
    p = Path(checkpoint_output) / f'{prefix}.{{arms}}.fasta' 
    glob_arms = glob_wildcards(p).arms
    res = expand(paths.telomere_variants.report, **wildcards, arms=glob_arms) ## create a list of inputs 
    return res


rule telomeric_kmers_agg:
    input: agg_telomeric_kmers_report
    output: paths.telomere_variants.final_report
    run:
        import pandas as pd
        from pathlib import Path
        
        def get_arm(f: str):
            f = Path(f)
            return f.name.split('.')[-3]

        df = pd.concat(
            [pd.read_csv(f,sep='\t').assign(arms=get_arm(f)) for f in input],
            ignore_index=True
        )
        df.to_csv(output[0], sep='\t', index=False)
        
        