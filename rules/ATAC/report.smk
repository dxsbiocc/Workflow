rule statsInfo:
    input:
        summary = expand("logs/bowtie2_{sample}.summary", sample=SAMPLES),
        spikein = expand("logs/bowtie2_{sample}_spikein.summary", sample=SAMPLES),
        metric = expand("dedup/{sample}/{sample}.metrics.txt", sample=SAMPLES),
        bam = expand("dedup/{pair}/{pair}.filtered.bam", pair=PAIRS),
        peak = expand("macs2/narrow/{pair}_peaks.narrowPeak", pair=PAIRS)
    output:
        "report/stats.csv"
    log:
        "logs/report_stats.log"
    params:
        sample_list = SAMPLES
    script:
        get_script("get_stats.py")