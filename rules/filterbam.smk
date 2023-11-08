
rule name_sorted_bam:
    input: 
        bam = paths.mapping.mmi_bam,
        idx = paths.mapping.mmi_bam_idx
    output: paths.mapping.name_sorted_bam
    conda: 'envs/minimap2.yml'
    shell:
        "samtools sort -n -@ 5 -m 4G -o {output} {input.bam} "


rule filterbam:
    input: 
        bam = paths.mapping.name_sorted_bam,
        # bam = paths.mapping.mmi_bam,
        # idx = paths.mapping.mmi_bam_idx,
        ncrf = paths.NoiseCancellingRepeatFinder.agg,
        tsv = lookup_value('anchor_tsv', mapping_df,'refgenome_id')
    output: 
        paths.filterbam.filtered_name_sorted
    shell:
        "python scripts/filter_bam_2023.py {input.bam} {input.ncrf} {input.tsv} {output}"
    
    # output: temp(paths.filterbam.sam)
    # shell:
    #     'python scripts/bam_filter_telomere.py {input.ncrf} {input.bed} {input.bam} {output}'

rule anchor_samtools_sort:
    input: paths.filterbam.filtered_name_sorted
    output: paths.filterbam.bam
    conda: 'envs/minimap2.yml'
    shell:
        'samtools sort {input} -o {output} '

rule anchor_samtools_index:
    input: paths.filterbam.bam
    output: paths.filterbam.bam_idx
    conda: 'envs/minimap2.yml'
    shell:
        'samtools index {input}  '

rule anchor_minimap2_stats:
    input: 
        bam = paths.filterbam.bam,
        idx = paths.filterbam.bam_idx
    output: paths.filterbam.stats
    params:
        tol = 50000,
    threads: 5
    # wrapper: 'file:wrappers/bam_stats'
    shell:
        'python scripts/bam_stats.v2.py '
        '-t {threads} '
        '--stdout '
        '-c '
        '{input.bam} > {output}'