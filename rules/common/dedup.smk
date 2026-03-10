if DEDUP == "markduplicates":
    rule markduplicates:
        input:
            bams = opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam"),
        output:
            bam = opj(OUTDIR, "dedup/{sample}/{sample}.rmdup.bam"),
            metrics = opj(OUTDIR, "dedup/{sample}/{sample}.metrics.txt"),
        log:
            opj(OUTDIR, "logs/dedup/markduplicates_{sample}.log"),
        params:
            extra = config['parameters']['picard']['markduplicates'],
        resources:
            mem_mb = config['parameters']['picard']['mem'],
        wrapper:
            get_wrapper('picard', 'markduplicates')
elif DEDUP == "sambamba":
    rule sambamba:
        input:
            opj(OUTDIR, "mapped/{sample}/{sample}.bam")
        output:
            opj(OUTDIR, "dedup/{sample}/{sample}.rmdup.bam")
        params:
            extra = ""  # optional parameters
        log:
            opj(OUTDIR, "logs/dedup/sambamba_markdup_{sample}.log")
        threads: 8
        shell:
            'sambamba view -f bam -F "not duplicate" {input}  | sambamba sort {params.extra} -o {output} /dev/stdin'
else:
    raise ValueError(f'the rule: {DEDUP} not support!')


