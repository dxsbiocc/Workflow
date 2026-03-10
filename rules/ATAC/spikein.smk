rule spikein:
    input:
        reads = expand(opj(OUTDIR, "trimmed/{{sample}}/{{sample}}.clean.{run}.fq.gz"), run=RUN),
        index = multiext(
            config['data']["ecoli"],
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        opj(OUTDIR, "mapped/{sample}/{sample}.spikein.bam"),
    log:
        opj(OUTDIR, "logs/mapped/bowtie2_{sample}_spikein.summary"),
    params:
        extra = config['parameters']['bowtie2']['extra'],  # optional parameters
        sort = config['parameters']['bowtie2']['sort']
    threads: config['parameters']['bowtie2']['threads']  # Use at least two threads
    wrapper:
        get_wrapper("bowtie2", "align")



# index bam
rule rmdup_index:
    input:
        opj(OUTDIR, "dedup/{sample}/{sample}.rmdup.bam"),
    output:
        opj(OUTDIR, "dedup/{sample}/{sample}.rmdup.bam.bai"),
    log:
        opj(OUTDIR, "logs/samtools/samtools_index_{sample}.log"),
    params:
        extra = "",  # optional params string
    wrapper:
        get_wrapper('samtools', 'index')
