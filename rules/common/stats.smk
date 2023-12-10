# stats for mapped bam
rule mapped_stats:
    input:
        bam = opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam")
    output:
        opj(OUTDIR, "report/mapped/stats/{sample}.stats")
    log:
        opj(OUTDIR, "logs/samtools/mapped_stats_{sample}.log"),
    wrapper:
        get_wrapper('samtools', 'stats')

rule mapped_plotBamStats:
    input:
        stats = rules.stats.output,
        gc = REF_GC
    output:
        directory(opj(OUTDIR, "report/mapped/plot/{sample}"))
    log:
        opj(OUTDIR, "logs/samtools/mapped_stats_plot_{sample}.log")
    wrapper:
        get_wrapper('samtools', 'plot-bamstats')
        

rule mapped_idxstats:
    input:
        bam = opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam")
    output:
        opj(OUTDIR, "report/mapped/stats/{sample}.idxstats")
    log:
        opj(OUTDIR, "logs/samtools/mapped_idxstats_{sample}.log"),
    params:
        extra = "",  # optional params string
    wrapper:
        get_wrapper('samtools', 'idxstats')


rule mapped_flagstat:
    input:
        bam = opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam")
    output:
        opj(OUTDIR, "report/mapped/stats/{sample}.flagstats")
    threads: 4
    log:
        opj(OUTDIR, "logs/samtools/mapped_flagstat_{sample}.log"),
    params:
        extra = "",  # optional params string
    wrapper:
        get_wrapper('samtools', 'flagstat')

# stats for dedupped bam
rule dedup_stats:
    input:
        bam = opj(OUTDIR, "dedup/{sample}/{sample}.rmdup.bam")
    output:
        opj(OUTDIR, "report/dedup/stats/{sample}.stats")
    log:
        opj(OUTDIR, "logs/samtools/dedup_stats_{sample}.log"),
    wrapper:
        get_wrapper('samtools', 'stats')

rule dedup_plotBamStats:
    input:
        stats = rules.stats.output,
        gc = REF_GC
    output:
        directory(opj(OUTDIR, "report/dedup/plot/{sample}"))
    log:
        opj(OUTDIR, "logs/samtools/dedup_stats_plot_{sample}.log")
    wrapper:
        get_wrapper('samtools', 'plot-bamstats')
        

rule dedup_idxstats:
    input:
        bam = opj(OUTDIR, "dedup/{sample}/{sample}.rmdup.bam")
    output:
        opj(OUTDIR, "report/dedup/stats/{sample}.idxstats")
    log:
        opj(OUTDIR, "logs/samtools/dedup_idxstats_{sample}.log"),
    params:
        extra = "",  # optional params string
    wrapper:
        get_wrapper('samtools', 'idxstats')


rule dedup_flagstat:
    input:
        bam = opj(OUTDIR, "dedup/{sample}/{sample}.rmdup.bam")
    output:
        opj(OUTDIR, "report/dedup/stats/{sample}.flagstats")
    threads: 4
    log:
        opj(OUTDIR, "logs/samtools/dedup_flagstat_{sample}.log"),
    params:
        extra = "",  # optional params string
    wrapper:
        get_wrapper('samtools', 'flagstat')