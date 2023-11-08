



rule ncrf:
    input: paths.basecalls.fasta
    output: paths.NoiseCancellingRepeatFinder.ncrf
    params:
        opts = get_opts(config['NCRF']['NCRF'], flag='--',fill_gap='=', skip_opts=['motif', 'exec_path']),
        opts_sort = get_opts(config['NCRF']['ncrf_sort']),
        motif = config['NCRF']['NCRF']['motif'],
        exec_path = config['NCRF']['NCRF']['exec_path'],
    conda: 'envs/ncrf.yml'
    shell:
        'cat {input} | {params.exec_path}/NCRF {params.motif} {params.opts}  > {output}'


rule ncrf_summary:
    input: paths.NoiseCancellingRepeatFinder.ncrf
    output: paths.NoiseCancellingRepeatFinder.summary
    params: 
        opts = '',
        exec_path = config['NCRF']['NCRF']['exec_path'],
    conda: 'envs/ncrf.yml'
    shell:
        "python2 --version && "
        'python2 {params.exec_path}/ncrf_cat.py {input} | python2 {params.exec_path}/ncrf_summary.py > {output}'


rule ncrf_summary_agg:
    input: paths.NoiseCancellingRepeatFinder.summary
    output: paths.NoiseCancellingRepeatFinder.agg
    shell:
        'python scripts/ncrf.summary.stats.py {input} {output}'

rule ncrf_reads_bed:
    input: paths.NoiseCancellingRepeatFinder.summary
    output: paths.NoiseCancellingRepeatFinder.bed
    shell:
        'python scripts/ncrf2bed.py {input} {output}'

rule telomeric_varients:
    input: 
        fasta = paths.telomeric.fasta,
        bed = paths.NoiseCancellingRepeatFinder.bed
    output: paths.NoiseCancellingRepeatFinder.variants
    shell:
        'python scripts/telomeric_variant_repeats.py {input} {output}'


rule list_of_telomeric_reads:
    '''List of telomeric candiadates '''
    input: paths.NoiseCancellingRepeatFinder.summary
    output: paths.telomeric.lor
    shell:
        '''awk -F '\\t' 'NR>1{{print $3}}' {input} | uniq > {output}'''


rule extract_telomeric_reads:
    '''extract telomeric read to fasta'''
    input:
        lor = paths.telomeric.lor,
        fasta = paths.basecalls.fasta
    output: paths.telomeric.fasta
    conda: 'envs/seqkit.yml'
    shell:
        'seqkit grep -f {input.lor} {input.fasta} > {output}'

rule extract_telomeres_only:
    '''extract only the telomere seq from reads'''
    input: 
        fa = paths.telomeric.fasta,
        bed = paths.NoiseCancellingRepeatFinder.bed
    output: paths.telomeric.telo_only
    conda: 'envs/bedtools.yml'
    shell:
        'bedtools getfasta -fi {input.fa} -fo {output} -bed {input.bed} '


##############################################################
##                alternative telomere variants             ##
##############################################################

def load_variants():
    fn = config['ALT_VARIANTS']['alt_motifs']
    motifs = []
    with open(fn) as fh:
        for motif in fh:
            motifs.append(motif.strip())
    return ' '.join(motifs)

rule alt_variants_ncrf:
    input: paths.telomeric.fasta
    output: paths.alt_variants.ncrf
    params:
        opts = get_opts(config['ALT_VARIANTS']['NCRF'], flag='--',fill_gap='=', skip_opts=['motif', 'exec_path']),
        opts_sort = get_opts(config['ALT_VARIANTS']['ncrf_sort']),
        motif = load_variants(), 
        exec_path = config['NCRF']['NCRF']['exec_path'],
    conda: 'envs/ncrf.yml'
    shell:
        'cat {input} | {params.exec_path}/NCRF {params.motif} {params.opts}  > {output}'


rule alt_variants_ncrf_summary:
    input: paths.alt_variants.ncrf
    output: paths.alt_variants.summary
    params: 
        opts = '',
        exec_path = config['NCRF']['NCRF']['exec_path'],
    conda:  'envs/ncrf.yml'
    shell:
        "python --version && "
        'python {params.exec_path}/ncrf_cat.py {input} | python {params.exec_path}/ncrf_summary.py > {output}'

