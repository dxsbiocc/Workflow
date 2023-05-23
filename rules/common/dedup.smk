rule markduplicates:
    input:
        bams = "mapped/{sample}/{sample}.sorted.bam",
    output:
        bam = "dedup/{sample}/{sample}.rmdup.bam",
        metrics = "dedup/{sample}/{sample}.metrics.txt",
    log:
        "logs/dedup/markduplicates_{sample}.log",
    params:
        extra = config['parameters']['picard']['markduplicates'],
    resources:
        mem_mb = 1024,
    wrapper:
        get_wrapper('picard', 'markduplicates')

rule sambamba:
    input:
        "mapped/{sample}/{sample}.bam"
    output:
        "dedup/{sample}/{sample}.rmdup.bam"
    params:
        extra = ""  # optional parameters
    log:
        "logs/dedup/sambamba_markdup_{sample}.log"
    threads: 8
    shell:
        'sambamba view -f bam -F "not duplicate" {input}  | sambamba sort {params.extra} -o {output} /dev/stdin'
