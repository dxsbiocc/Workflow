if not os.path.exists(INDEX):
    rule simpleaf_index:
        input:
            fasta = REFERENCE,
            gtf = GTF
        output:
            directory(INDEX),
        threads: 20
        params:
            extra = "-r 90",
        log:
            "logs/quant/simpleaf_index.log"
        wrapper:
            get_wrapper("simpleaf", "index")

rule simpleaf_quant:
    input:
        fq1 = lambda wildcards: DATA.loc[wildcards.sample, 'fastq1'].to_list(),
        fq2 = lambda wildcards: DATA.loc[wildcards.sample, 'fastq2'].to_list(),
        index = INDEX,
    output:
        directory(opj(OUTDIR, "simpleaf/{sample}")),
    log:
        "logs/quant/simpleaf_{sample}.log",
    params:
        extra = "-c 10xv2 -u ",
        resolution = "cr-like"  # possible values: cr-like, cr-like-em, parsimony, parsimony-em, parsimony-gene, parsimony-gene-em
    threads: 8
    wrapper:
        get_wrapper("simpleaf", "quant")