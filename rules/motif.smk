rule get_peak:
    input:
        "macs2/narrow/{sample}_peaks.narrowPeak"
    output:
        "macs2/narrow/{sample}_homer_peaks.txt"
    shell:
        """awk '{{print $4"\t"$1"\t"$2"\t"$3"\t""+"}}' {input} > {output}"""
    

rule find_motifs:
    input:
        peak = "macs2/narrow/{sample}_homer_peaks.txt",
        genome = config['data']['ref']
    output:
        protected("motifs/{sample}/")
    params:
        extra = "-size 200"
    threads: 2
    log:
        "logs/findMotifs_{sample}.log"
    wrapper:
        get_wrapper("findMotifsGenome")