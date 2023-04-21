rule markduplicates:
    input:
        bams = "mapped/{sample}.sorted.bam",
    output:
        bam = "dedup/{sample}.rmdup.bam",
        metrics = "dedup/{sample}.metrics.txt",
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
        "mapped/{sample}.bam"
    output:
        "dedup/{sample}.rmdup.bam"
    params:
        extra = ""  # optional parameters
    log:
        "logs/dedup/sambamba_markdup_{sample}.log"
    threads: 8
    shell:
        'sambamba view -f bam -F "not duplicate" {input}  | sambamba sort {params.extra} -o {output} /dev/stdin'
