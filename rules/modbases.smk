
print(basecalls_df)
use rule minimap2_mapping as modbase_mapping with:
    input:
        ref = paths.refgenome.fasta,
        idx = paths.refgenome.mmi_idx,
        fq = lookup_value('mod_fastq', basecalls_df)
    output: paths.modbases.bam

rule merge_bam:
    input: expand_rows_with_lookup(paths.modbases.bam, mapping_df, is_aggregated=True )
    output: paths.modbases.merged_bam
    shell:
        "samtools merge -o {output} {input}"


use rule minimap2_index as modbase_mapped_index with :
    input: paths.modbases.merged_bam
    output: paths.modbases.bam_idx

rule extact_mods:
    input:
        ref = paths.refgenome.fasta,
        bam = paths.modbases.merged_bam,
        bam_idx = paths.modbases.bam_idx,
        # bed = paths.refgenome.bed,
    output:
        paths.modbases.bed
    params:
        path_exec = config['MODKIT']['path_exec'],
        opt = get_opts(config['MODKIT'], flag='mix', skip_opts=['path_exec'])
    shell:
        "{params.path_exec} pileup {input.bam} {output} --ref {input.ref} {params.opt}"
        

rule plot_mods:
    input:
        methyl = paths.modbases.bed,
        bed = paths.refgenome.bed,
    output:
        pdf = paths.modbases.pdf,
        csv = paths.modbases.csv
    params:
        min_cov = 10
    shell:
        "python scripts/plot_methlyation_from_modkit.py --min-read-cov {params.min_cov} {input.bed} {input.methyl} {output.pdf} "