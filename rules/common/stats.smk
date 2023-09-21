rule stats:
    input:
        bam = "dedup/{sample}/{sample}.rmdup.bam"
    output:
        "report/stats/{sample}.stats"
    log:
        "logs/samtools/stats_{sample}.log",
    wrapper:
        get_wrapper('samtools', 'stats')

rule plotBamStats:
    input:
        stats = rules.stats.output,
        gc = REF_GC
    output:
        directory("report/plot/{sample}")
    log:
        "logs/samtools/stats_plot_{sample}.log"
    wrapper:
        get_wrapper('samtools', 'plot-bamstats')
        

rule idxstats:
    input:
        bam = "dedup/{sample}/{sample}.rmdup.bam"
    output:
        "report/stats/{sample}.idxstats"
    log:
        "logs/samtools/idxstats_{sample}.log",
    params:
        extra = "",  # optional params string
    wrapper:
        get_wrapper('samtools', 'idxstats')


rule flagstat:
    input:
        bam = "dedup/{sample}/{sample}.rmdup.bam"
    output:
        "report/stats/{sample}.flagstats"
    threads: 4
    log:
        "logs/samtools/flagstat_{sample}.log",
    params:
        extra = "",  # optional params string
    wrapper:
        get_wrapper('samtools', 'flagstat')