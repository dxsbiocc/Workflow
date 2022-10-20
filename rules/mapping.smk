# align to reference genome
rule bowtie2:
    input:
        reads = [
            "trimmed/{sample}.clean.R1.fq.gz", 
            "trimmed/{sample}.clean.R2.fq.gz"
        ],
        index = multiext(
            config['data']["index"],
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        "bowtie2/{sample}.sort.bam",
    log:
        "logs/bowtie2_{sample}.log",
    params:
        extra = config['parameters']['bowtie2']['extra'],  # optional parameters
        sort = config['parameters']['bowtie2']['sort']
    threads: 8  # Use at least two threads
    wrapper:
        get_wrapper("bowtie2", "align")

# remove dupliactes
rule mark_duplicates:
    input:
        bams = protected("bowtie2/{sample}.sort.bam"),
    # optional to specify a list of BAMs; this has the same effect
    # of marking duplicates on separate read groups for a sample
    # and then merging
    output:
        bam = "dedup/{sample}.cons.bam",
        metrics = "dedup/{sample}.metrics.txt",
    log:
        "logs/picard_dedup_{sample}.log",
    params:
        extra = config['parameters']['picard']['MarkDuplicates'],
    resources:
        mem_mb = 1024,
    wrapper:
        get_wrapper('markduplicates')

# index bam
rule cons_index:
    input:
        "dedup/{sample}.cons.bam",
    output:
        "dedup/{sample}.cons.bam.bai",
    log:
        "logs/samtools_index_{sample}.log",
    params:
        extra = "",  # optional params string
    wrapper:
        get_wrapper('samtools', 'index')
