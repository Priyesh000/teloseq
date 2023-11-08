
## combind bam sam species data
## bed region (ex/inculding telomere)
## call variants -> vcf
## phase -> bam

## New steps:
    ## Filter for full length
    ## Checkup point for medaka
    


def _get_region_by_coverage(file_path, flag=None, shift_one_base=False ):

    df = pd.read_csv(file_path, sep='\t')
    df = df[df.passed] 
    ##BedFmt.from_region
    list_of_region = []
    for region in df.region.to_list():
        b = BedFmt.from_region(region)
        if shift_one_base:
            b.shift_one_base ## from zero base index to one base
        if flag:
            region = f'{flag} ' + b.to_region()
        list_of_region.append(region) 
    return list_of_region        



rule vcf_compress:
    input: paths.variants.hapcut2_block_vcf,
    output: paths.variants.hapcut2_block_vcf_bgzip
    conda: 'envs/medaka.yml'
    shell:
        # runs bgzip on the input and tabix indexes the result
        ##bgzip >[file] && tabix -fp vcf [file]
        'cat {input} | vcfstreamsort | bgziptabix {output}'


rule whatshap_haplotag:
    input:
        vcf = paths.variants.hapcut2_block_vcf_bgzip,
        bam = paths.variants.bam,
        ref = paths.refgenome.fasta
    output: 
        bam = paths.variants.phased_bam_hapcut,
        # hap = paths.variants.haplotag 
    params:
        opts = get_opts(config['WHATSHAP']['whatshap_tag'], flag='mix'),
        ## if regions not given then if fail on chromosomes without phasing
        regions = lambda wildcards: _get_region_by_coverage(paths.variants.coverage.format(**wildcards), flag='--regions', shift_one_base=True)
    conda: 'envs/medaka.yml'
    shell:
        'whatshap haplotag {params.opts} -o {output.bam} '
        # '--output-haplotag-list {output.hap} '  
        '{params.regions} ' 
        '--reference {input.ref}  {input.vcf} {input.bam} ;  '
        'sleep 30; '
        'samtools index {output}'

