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
        "logs/{sample}_bowtie2.summary",
    params:
        extra = config['parameters']['bowtie2']['extra'],  # optional parameters
        sort = config['parameters']['bowtie2']['sort']
    threads: 8  # Use at least two threads
    wrapper:
        get_wrapper("bowtie2", "align")

rule spikein:
    input:
        reads = [
            "trimmed/{sample}.clean.R1.fq.gz", 
            "trimmed/{sample}.clean.R2.fq.gz"
        ],
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
        "bowtie2/{sample}.spikein.bam",
    log:
        "logs/{sample}_bowtie2_spikein.summary",
    params:
        extra = "--end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 -I 10 -X 700",  # optional parameters
        sort = "none"
    threads: 8  # Use at least two threads
    wrapper:
        get_wrapper("bowtie2", "align")

rule stats:
    input:
        bam = rules.bowtie2.output
    output:
        "bowtie2/stats/{sample}.stats"
    log:
        "logs/samtools_stats_{sample}.log",
    wrapper:
        get_wrapper('samtools', 'stats')

rule plotBamStats:
    input:
        stats = rules.stats.output,
        gc = GC
    output:
        directory("bowtie2/plot/{sample}")
    log:
        "logs/samtools_stats_plot_{sample}.log"
    wrapper:
        get_wrapper('samtools', 'plot-bamstats')

rule idxstats:
    input:
        bam = rules.bowtie2.output
    output:
        "bowtie2/stats/{sample}.idxstats"
    log:
        "logs/samtools_idxstats_{sample}.log",
    params:
        extra = "",  # optional params string
    wrapper:
        get_wrapper('samtools', 'idxstats')

rule flagstat:
    input:
        bam = rules.bowtie2.output
    output:
        "bowtie2/stats/{sample}.flagstats"
    threads: 4
    log:
        "logs/samtools_flagstat_{sample}.log",
    params:
        extra = "",  # optional params string
    wrapper:
        get_wrapper('samtools', 'flagstat')

# remove dupliactes
rule mark_duplicates:
    input:
        bams = "bowtie2/{sample}.sort.bam",
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
        get_wrapper('picard', 'markduplicates')

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
