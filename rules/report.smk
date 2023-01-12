rule statsInfo:
    input:
        summary = expand("logs/{sample}_bowtie2.summary", sample=SAMPLES),
        spikein = expand("logs/{sample}_bowtie2_spikein.summary", sample=SAMPLES),
        metric = expand("dedup/{sample}.metrics.txt", sample=SAMPLES),
        bam = expand("dedup/{sample}.filtered.bam", sample=SAMPLES),
        peak = expand("macs2/narrow/{pair}_peaks.narrowPeak", pair=PAIRS)
    output:
        "report/stats.csv"
    params:
        sample_list = SAMPLES
    script:
        get_script("get_stats.py")