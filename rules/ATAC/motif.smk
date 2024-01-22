rule get_peak:
    input:
        get_peakfile
    output:
        opj(OUTDIR, "motifs/anno/{pair}_homer_peaks.txt")
    shell:
        """awk '{{print $4"\t"$1"\t"$2"\t"$3"\t""+"}}' {input} > {output}"""

rule find_motifs:
    input:
        peak = rules.get_peak.output,
        genome = config['data']['ref']
    output:
        directory(opj(OUTDIR, "motifs/{pair}/"))
    params:
        extra = "-size 200"
    threads: 2
    log:
        opj(OUTDIR, "logs/motifs/findMotifs_{pair}.log")
    wrapper:
        get_wrapper('homer', "findMotifsGenome")