rule statsInfo:
    input:
        summary = expand(opj(OUTDIR, "logs/bowtie2_{sample}.summary"), sample=SAMPLES),
        spikein = expand(opj(OUTDIR, "logs/bowtie2_{sample}_spikein.summary"), sample=SAMPLES),
        metric = expand(opj(OUTDIR, "dedup/{sample}/{sample}.metrics.txt"), sample=SAMPLES),
        bam = expand(opj(OUTDIR, "dedup/{pair}/{pair}.filtered.bam"), pair=PAIRS),
        peak = expand(opj(OUTDIR, "macs2/narrow/{pair}_peaks.narrowPeak"), pair=PAIRS)
    output:
        opj(OUTDIR, "report/stats.csv")
    log:
        opj(OUTDIR, "logs/report_stats.log")
    params:
        sample_list = SAMPLES
    script:
        get_script("get_stats.py")