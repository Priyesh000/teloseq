rule TVR_signal:
    input: 
        _ = paths.merged_bams.bam_idx,
        bam = paths.merged_bams.bam,
        bed = paths.refgenome.tsv
    output: 
        png = paths.TVR.png,
        # svg = paths.TVR.svg, 
        csv = paths.TVR.csv,
        # raw = paths.TVR.raw
    params: 
        min_telomere_length = 3000,
        outdir = lambda w, output: Path(output.png).parent,
        prefix = lambda w, output: Path(output.png).stem.replace(".tvr_signal.norm", "")
    shell:
        "echo {params.prefix} && "
        "python scripts/TVR_signal.py --bamfile {input.bam} "
        "--anchor_bed {input.bed} "
        "--output_dir {params.outdir} "
        "--prefix {params.prefix} "
        "--min_seqeunce_length {params.min_telomere_length} "


rule copy2report:
    input: 
        png = paths.TVR.png,
        csv = paths.TVR.csv
    output: 
        png = paths.reports.png,
        csv = paths.reports.csv
    shell:
        "cp {input.png} {output.png} && cp {input.csv} {output.csv}"
