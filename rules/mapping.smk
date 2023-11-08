
rule minimap2_ref_indx:
    input: paths.refgenome.fasta
    output: paths.refgenome.mmi_idx
    conda: 'envs/minimap2.yml'
    params:
        opts = get_opts(config['MINIMAP2']['index'], flag='-', skip_opts=['t'])    
    shell: 
        'minimap2 {params.opts} -d {output} {input}'

rule minimap2_mapping:
    '''Align with minimap2'''
    input: 
        ref = paths.refgenome.fasta,
        idx = paths.refgenome.mmi_idx,
        fq = paths.basecalls.fasta
    output: paths.mapping.mmi_bam
    params:
        opts = get_opts(config['MINIMAP2']['mapping'], flag='-', skip_opts=['t'])    
    threads: config['MINIMAP2']['mapping']['t']
    conda: 'envs/minimap2.yml'
    shell:
        'minimap2 -t {threads} {params.opts} {input.idx} {input.fq} | samtools view -h --reference {input.ref} | samtools sort -o {output} '

rule minimap2_index:
    input: paths.mapping.mmi_bam
    output: paths.mapping.mmi_bam_idx
    conda: 'envs/samtools.yml'
    shell:
        'samtools index {input}'

rule minimap2_stats:
    input: 
        bam = paths.mapping.mmi_bam,
        _ = paths.mapping.mmi_bam_idx
    output: paths.mapping.mmi_stat
    params:
        tol = 25000
    threads: config['MINIMAP2']['mapping']['t']
    # wrapper: 'file:wrappers/bam_stats'
    shell:
        'python scripts/bam_stats.v2.py '
        '-t {threads} '
        '--stdout '
        '-c '
        '{input.bam} > {output}'
    

# rule winnowmap_stats:
#     input: paths.mapping.wmap_bam
#     output: paths.mapping.wmap_stat
#     #conda: 'env/general.yml'
#     params:
#         tol = 50000,
#     threads: config['WINNOWMAP']['mapping']['t']
#     wrapper: 'file:wrappers/bam_stats'
#     # shell:
#     #     'python scripts/bam_stats.py  -t {threads} --tol 50000 -s  {input} > {output}'


rule create_bed_file:
    input:  paths.refgenome.tsv
    output: paths.refgenome.bed
    shell:
        "python scripts/build_chrom_ends_bedfile.py {input} {output}"
        

rule mosdepth_minimap_ends_coverage:
    input: 
        bam = paths.mapping.mmi_bam,
        _   = paths.mapping.mmi_bam_idx,
        bed = paths.refgenome.bed
    output: paths.mapping.mmi_cov
    params:
        prefix = lambda x,output: output[0].replace('.mosdepth.global.dist.txt', ''),
        opts = get_opts(config['MOSDEPTH'], flag='mix', skip_opts=['t','threads', 'prefix', ]),
    conda: 'envs/mosdepth.yml'
    threads: config['MOSDEPTH']['threads']
    shell:
        'mosdepth {params.opts} -t {threads} --by {input.bed} {params.prefix} {input.bam}'


# rule mosdepth_winnowmap_windows:
#     input: 
#         bam = paths.mapping.wmap_bam,
#         _   = paths.mapping.wmap_bam_idx,
#     output: paths.coverage_check.wmap_window
#     params:
#         prefix = lambda x,output: output[0].replace('.regions.bed.gz', ''),
#         opts = get_opts(config['MOSDEPTH'], flag='mix', skip_opts=['t','threads', 'prefix', ]),
#         by = 1000
#     conda: 'envs/mosdepth.yml'
#     threads: config['MOSDEPTH']['threads']
#     shell:
#         'mosdepth {params.opts} -t {threads} --by {params.by} {params.prefix} {input.bam}'

# rule mosdepth_minimap_windows:
#     input: 
#         bam = paths.mapping.mmi_bam,
#         _   = paths.mapping.mmi_bam_idx,
#     output: paths.coverage_check.mmi_window
#     params:
#         prefix = lambda x,output: output[0].replace('.regions.bed.gz', ''),
#         opts = get_opts(config['MOSDEPTH'], flag='mix', skip_opts=['t','threads', 'prefix', ]),
#         by = 50000
#     conda: 'envs/mosdepth.yml'
#     threads: config['MOSDEPTH']['threads']
#     shell:
#         'mosdepth {params.opts} -t {threads} --by {params.by} {params.prefix} {input.bam}'



# rule mosdepth_lra:
#     input: 
#         bam = paths.mapping.lra_bam,
#         _   = paths.mapping.lra_bam_idx,
#         bed = paths.refgenome.restricted
#     output: paths.mapping.lra_cov
#     params:
#         prefix = lambda x, output : output[0].replace('.mosdepth.global.dist.txt', ''),
#         opts = get_opts(config['MOSDEPTH'], flag='mix', skip_opts=['t','threads', 'prefix', ]),
#         # by = './scripts/restricted.chromosome_arms.20Kb.bed'
#     conda: 'envs/mosdepth.yml'
#     threads: config['MOSDEPTH']['threads']
#     shell:
#         'mosdepth {params.opts} -t {threads} --by {input.bed} {params.prefix} {input.bam}'

# rule window_mosdepth_lra:
#     input: 
#         bam = paths.mapping.lra_bam,
#         _   = paths.mapping.lra_bam_idx,
#     output: paths.coverage_check.lra_window
#     params:
#         prefix = lambda x, output : output[0].replace('.regions.bed.gz', ''),
#         opts = get_opts(config['MOSDEPTH'], flag='mix', skip_opts=['t','threads', 'prefix' ]),
#     conda: 'envs/mosdepth.yml'
#     threads: config['MOSDEPTH']['threads']
#     shell:
#         'mosdepth {params.opts} -t {threads} --by 100000 {params.prefix} {input.bam}'

# rule window_mosdepth_lra_add_arms:
#     input: paths.coverage_check.lra_window
#     output: paths.coverage_check.lra_window_arms
#     conda: 'envs/bedtools.yml'
#     shell:
#         'bedtools intersect -a {input} -b scripts/arms_transition_coord.bed -bed -wao | cut -f1-4,8 > {output}'


# rule blastn_db:
#     input: 
#         paths.refgenome.rs_ncrf_masked
#         # paths.refgenome.fasta
#     output: paths.refgenome.blast_idx
#     conda: 'envs/blastn.yml'    
#     shell:
#         'export BLASTDB_LMDB_MAP_SIZE=1000000; '
#         'makeblastdb -in {input} -dbtype nucl'

# rule blastn:
#     input: 
#         # ref = paths.refgenome.fasta,
#         ref = paths.refgenome.rs_ncrf_masked,
#         _ = paths.refgenome.blast_idx,
#         tela = paths.telomeric.fasta
#     output: paths.blast.align
#     params:
#         opts = get_opts(config['BLASTn']['blastn'], flag='-', skip_opts=['num_threads'])
#     conda: 'envs/blastn.yml'
#     threads: config['BLASTn']['blastn']['num_threads']
#     shell:
#         'blastn -db {input.ref} -query {input.tela} {params.opts} -num_threads {threads} -out {output}'        

def get_subtelomere_span(wildcards):
    reference = config['FULL_LENGTH']['subtelomere-span']
    ref_id = wildcards.get('refgenome_id')
    print(ref_id)
    return reference.get(ref_id, 500)


# rule full_length:
#     input: 
#         bam = paths.mapping.mmi_bam,
#         _ =  paths.mapping.mmi_bam_idx
#     output: paths.mapping.mmi_fl
#     params:
#         telomere_bed = lambda wildcards: get_column_value('telomere_bed_path', ref_df, wildcards.get('refgenome_id')),
#         subtelo_span = lambda wildcards: get_subtelomere_span(wildcards),
#         opts = get_opts(config['FULL_LENGTH'], flag='--', skip_opts=['subtelomere-span'])
#     shell:
#         'python scripts/full_lengeth_mapped_reads.py '
#         '--bed {params.telomere_bed} '
#         '--bam {input.bam} '
#         '--output {output}.tmp '
#         '--subtelomere-span {params.subtelo_span} '
#         '{params.opts} '
#         '&& samtools sort -o {output} {output}.tmp '
#         '&& samtools index {output} '
#         '&& rm -fr {output}.tmp '


# rule full_length_stats:
#     input: paths.mapping.mmi_fl
#     output: paths.mapping.mmi_fl_stats
#     params:
#         tol = 50000,
#     threads: config['MINIMAP2']['mapping']['t']
#     wrapper: 'file:wrappers/bam_stats'


rule merge_bams_files:
    input: expand_rows_with_lookup(paths.filterbam.bam, mapping_df, is_aggreated=True)
    output: paths.merged_bams.bam
    conda: 'envs/samtools.yml'
    threads: config['MINIMAP2']['mapping']['t']
    shell:
        'samtools merge -f --threads {threads}  {output} {input}'

# use rule minimap2_index as merged_bam_index with:
#     input: paths.merged_bams.bam
#     output: paths.merged_bams.bam_idx
rule merged_bam_index:
    input: paths.merged_bams.bam
    output: paths.merged_bams.bam_idx
    conda: 'envs/samtools.yml'
    shell:
        'samtools index {input}'