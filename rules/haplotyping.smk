rule merge_telomeres_fasta:
    input: expand_rows_with_lookup(paths.telomeric.fasta, basecalls_df, is_aggregated=True)
    output: paths.merged_telomeres.fasta
    shell:
        "cat {input} > {output}"

def get_ref():
    res = ref_df.loc['hg002_paternal', 'refgenome_path']
    tsv = ref_df.loc['hg002_paternal', 'anchor_tsv']
    path = Path(res)
    if not path.exists():
        raise ValueError(f'File does not exists {path}')
    return path, tsv

rule copy_ref:
    input: get_ref()[0]
    output: 
        ref = paths.wf_human_var.ref
    shell:
        'cp {input} {output}'


use rule copy_ref as copy_ref_tsv with:
    input: get_ref()[1] 
    output: paths.wf_human_var.ref_tsv


use rule minimap2_ref_indx as ref_paternal_index with:
    input: paths.wf_human_var.ref
    output: paths.wf_human_var.ref_mmi

use rule minimap2_mapping as map_paternal with:
    input: 
        # ref = lookup_value('refgenome_path', ref_df),
        ref = paths.wf_human_var.ref,
        idx = paths.wf_human_var.ref_mmi,
        fq = paths.merged_telomeres.fasta
    output: paths.wf_human_var.bam 

use rule minimap2_index as map_index with:
    input: paths.wf_human_var.bam
    output: paths.wf_human_var.bam_idx

rule samtools_faidx:
    input: paths.wf_human_var.ref
    output: paths.wf_human_var.ref_idx
    conda: 'envs/samtools.yml' 
    shell:
        'samtools faidx {input}'

rule download_cliarmodel:
    output: directory(paths.models.clair)
    params: inputfile = config['CLAIR3']['url']
    # log: to_log(paths.models.clair)
    wrapper:
        'file://wrappers/downloader/'

def get_wf_human_var_output_files():
    data = config["WF_HUMAN_VARIATION"]
    vaild_options = ['snp', 'sv', 'methyl', 'cnv', 'str']
    output_files = []
    for key, val in data.items():
        if val is True:
            filepath = getattr(paths.wf_human_var, key, None)
            if filepath:
                output_files.append(filepath)
    # print(output_files)
    return output_files

rule telomere_regions:
    input: paths.wf_human_var.ref_tsv
    output: paths.wf_human_var.ref_bed
    shell:
        "awk -F'\t' 'NR>1{{print $1\"\t\"$2\"\t\"$3}}' {input} > {output}"

rule run_clair3:
    input: 
        bam = paths.wf_human_var.bam,
        _ = paths.wf_human_var.bam_idx,
        ref = paths.wf_human_var.ref,
        _ref_idx = paths.wf_human_var.ref_idx,
        bed_region = paths.wf_human_var.ref_bed,
        # model_path = paths.models.clair,
    output: 
        paths.wf_human_var.snp,
        paths.wf_human_var.phased_snp,
        paths.wf_human_var.phased_bam
    params: 
        platform = "ont",
        model_path = paths.models.clair,
        output = lambda w, output: (Path(output[0]).resolve().parent) 
    conda: 'envs/clair3.yml'
    threads: 30
    shell:
        "rm -fr {params.output} && "
        "mkdir -p {params.output}/tmp/split_beds/ && "
        "run_clair3.sh --bam_fn={input.bam} "
        "--ref_fn={input.ref} "
        "--output={params.output} "
        "--platform={params.platform} "
        "--model_path={params.model_path} "
        "--bed_fn={input.bed_region} "
        "--threads={threads} --chunk_size=25000 "
        "--use_whatshap_for_final_output_haplotagging "
        
        # "--use_longphase_for_final_output_phasing && "
        # "--use_whatshap_for_final_output_phasing "
        # "&& "
        # "if [[ -e {params.output} ]]; then "
