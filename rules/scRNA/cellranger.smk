rule cellranger_index:
    input:
        fasta = REFERENCE,
        gtf = GTF
    output:
        directory(INDEX)
    threads: 8
    log:
        opj(OUTDIR, "logs/count/cellranger_index.log"),
    params:
        genome = GENOME,
        cellranger = CELLRANGER
    wrapper:
        get_wrapper("cellranger", "mkref")

rule cellranger_count:
    input:
        index = INDEX,
        fastqs = expand("data/{{sample}}.{lane}.fastq", lane=SAMPLES),
    output:
        directory(opj(OUTDIR, "cellranger/{sample}"))
    threads: 8
    log:
        opj(OUTDIR, "logs/count/cellranger_count.log"),
    params:
        cellranger = CELLRANGER,
    wrapper:
        get_wrapper("cellranger", "count")