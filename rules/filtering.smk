rule rm_chrM:
    input:
        "dedup/{sample}.cons.bam",
    output:
        "dedup/{sample}.rmChrM.bam"
    log:
        "logs/dedup_rmChrM_{sample}.log"
    run:
        shell("(samtools view -h {input} | grep -v chrM | samtools view -bS | samtools sort -@ 4 -o {output})")
        shell("(samtools index {output})")

rule shift:
    input:
         "dedup/{sample}.rmChrM.bam"
    output:
        "dedup/{sample}.shift.bam"
    log:
        "logs/dedup_{sample}_shift.log"
    threads:
        4
    params:
        # optional parameters
        extra = ""
    wrapper:
        get_wrapper("deeptools", "alignmentSieve")

rule shift_sort:
    input:
        "dedup/{sample}.shift.bam"
    output:
        "dedup/{sample}.shift.sort.bam"
    log:
        "log/dedup_{sample}_shift.sort.log",
    params:
        extra="-m 4G",
    threads: 8
    wrapper:
        get_wrapper('samtools', 'sort')

rule shift_index:
    input:
        "dedup/{sample}.shift.sort.bam",
    output:
        "dedup/{sample}.shift.sort.bam.bai",
    log:
        "logs/dedup_{sample}_shift_index.log",
    params:
        extra = "",  # optional params string
    wrapper:
        get_wrapper('samtools', 'index')