rule statsInfo:
    input:
        summary = expand("logs/{sample}_bowtie2.summary", sample=samples.index),
        spikein = expand("logs/{sample}_bowtie2_spikein.summary", sample=samples.index),
        metric = expand("dedup/{sample}.metrics.txt", sample=samples.index),
        bam = expand("dedup/{sample}.filtered.bam", sample=samples.index),
        peak = expand("macs2/narrow/{sample}_peaks.narrowPeak", sample=samples.index)
    output:
        "report/stats.csv"
    params:
        sample_list = samples.index.to_list()
    script:
        get_script("get_stats.py")