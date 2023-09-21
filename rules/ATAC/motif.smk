rule get_peak:
    input:
        "macs2/narrow/{pair}_peaks.narrowPeak"
    output:
        "macs2/narrow/{pair}_homer_peaks.txt"
    shell:
        """awk '{{print $4"\t"$1"\t"$2"\t"$3"\t""+"}}' {input} > {output}"""

rule find_motifs:
    input:
        peak = "macs2/narrow/{pair}_homer_peaks.txt",
        genome = config['data']['ref']
    output:
        directory("motifs/{pair}/")
    params:
        extra = "-size 200"
    threads: 2
    log:
        "logs/findMotifs_{pair}.log"
    wrapper:
        get_wrapper('homer', "findMotifsGenome")