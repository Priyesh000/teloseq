import pandas as pd
from pathlib import Path
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.4.0")

## wildcard_constraints:
#
##### load config and sample sheets #####

configfile: "config.yml"

##### load rules #####

include: "rules/common.smk"

basecalls_df, ref_df, mapping_df = create_input_dataframes()
paths  = create_path_accessor(config['FILE_LAYOUT'], Path(config['SAVEPATH']) )

# print(mapping_df.columns)
include: "rules/ingest.smk"
include: "rules/repeat_finder.smk"
include: "rules/mapping.smk"
# include: "rules/kmers.smk"
include: "rules/telomere_variants.smk"
include: "rules/stats.smk"
include: "rules/filterbam.smk"
include: "rules/reports.smk"
include: "rules/modbases.smk"
include: "rules/tvr.smk"
include: "rules/haplotyping.smk"

# constrain wildcards
wildcard_constraints:
    dataset_id = "(" + '|'.join(basecalls_df.dataset_id.unique()) + ")", ###'(HelaVST|HelaLT|HG002|U2OS|)',
    refgenome_id = "(" + '|'.join(ref_df.refgenome_id.unique()) + ")"
    # run_id = '([P|F][A-Z0-9]+)'


hap_df = mapping_df.copy()
index_name = hap_df.index.names
hap_df.loc[:, 'refgenome_id'] = 'hg002_paternal'
hap_df.set_index(index_name, inplace=True, drop=False)


##### target rules #####
rule all:
    """Full pipeline"""
    input: 
        expand_rows(paths.basecalls.fasta, basecalls_df),
        expand_rows(paths.basecalls.stats, basecalls_df),
        expand_rows(paths.basecalls.agg_stats, basecalls_df),
        # expand_rows(paths.basecalls.summary, basecalls_df),
        expand_rows(paths.NoiseCancellingRepeatFinder.summary, basecalls_df),
        expand_rows(paths.mapping.mmi_bam_idx, mapping_df),
        expand_rows(paths.mapping.mmi_stat, mapping_df),
        expand_rows(paths.mapping.mmi_cov, mapping_df),

        expand_rows(paths.NoiseCancellingRepeatFinder.variants, mapping_df),
        expand_rows(paths.NoiseCancellingRepeatFinder.agg, mapping_df),
        # expand_rows(paths.mapping.mmi_fl_stats, mapping_df),
        expand_rows(paths.stats.agg_stat, basecalls_df),
        expand_rows(paths.filterbam.stats, mapping_df),
        expand_rows(paths.merged_bams.bam, mapping_df),
        
        ## TODO: reports on merge bams
        expand_rows(paths.reports.stats_csv, mapping_df),
        expand_rows(paths.reports.report_merge, mapping_df),
        expand_rows(paths.reports.html, mapping_df),
        # expand_rows(paths.TVR.png, mapping_df),  
        # expand_rows(paths.TVR.csv, mapping_df),

rule modbases:
    """Modbases on filtered and merged bams"""
    input:
        expand_rows(paths.modbases.bed, mapping_df),
        expand_rows(paths.modbases.pdf, mapping_df),
        expand_rows(paths.modbases.csv, mapping_df)

rule TVR:
    """Telomere Variant Repeats on filtered bam and merged bams"""
    input:
        expand_rows(paths.reports.png, mapping_df),  
        expand_rows(paths.reports.csv, mapping_df)
        
rule hap:
    """Haplotype telomeric reads"""
    input:
        expand_rows(paths.wf_human_var.bam , hap_df),
        expand_rows(paths.wf_human_var.snp, hap_df),


# rule var:
#     """Identify telomeric variance"""
#     input: 
#         expand_rows(paths.telomeric.telo_only, basecalls_df),
#         expand_rows(paths.telomeric.fa_merged, mapping_df),
#         # expand_rows(paths.telomere_variants.bam_merged, mapping_df),
#         expand_rows(paths.telomere_variants.final_report, mapping_df)
