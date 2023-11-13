
rule report_per_run:
    input: 
        ncrf = paths.NoiseCancellingRepeatFinder.agg,
        bcalls = paths.basecalls.stats,
        bam = paths.filterbam.stats
    output: paths.reports.stats_csv
    # output: paths
    conda: "envs/jupyter.yml"
    log:
        notebook=paths.reports.report 
    notebook:
        '../report_template/20221130_telomere_report_template.py.ipynb'

rule report_to_html:
    input: paths.reports.report
    output: paths.reports.html
    params:
        dirpath = lambda w, output: Path(output[0]).parent,
        filename = lambda w, output: Path(output[0]).name,
    conda: "envs/jupyter.yml"
    shell:
        'jupyter nbconvert --to html --template lab --no-input {input} --output {params.filename} --output-dir {params.dirpath}'


rule merge_telomeres_stats:
    input: 
        expand_rows_with_lookup(paths.NoiseCancellingRepeatFinder.agg, basecalls_df, is_aggregated=True)
    output: paths.reports.ncrf
    run:
        import polars as pl 
        if len(input) == 1:
            pl.read_csv(input[0]).write_csv(output[0])
        else:
            pl.concat([pl.scan_csv(str(f)) for f in input]).collect().write_csv(output[0])

use rule merge_telomeres_stats as merge_basecall_stats with:
    input: 
        expand_rows_with_lookup(paths.basecalls.stats, basecalls_df, is_aggregated=True)
    output: paths.reports.basecall_stats

use rule merge_telomeres_stats as merge_bam_stats with:
    input: 
        expand_rows_with_lookup(paths.filterbam.stats, mapping_df, is_aggregated=True)
    output: paths.reports.bam_stats 

use rule report_per_run as report_per_sample with:
    input: 
        ncrf = paths.reports.ncrf,
        bcalls = paths.reports.basecall_stats,
        bam = paths.reports.bam_stats
    output: paths.reports.stats_csv_merged
    log:
        notebook=paths.reports.report_merge
