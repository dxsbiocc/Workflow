rule baserecalibrator:
    input:
        bam = "dedup/{sample}.rmdup.bam",
        ref = REFERENCE,
        dict = "resources/genome.dict",
        known = KNOWN,
    output:
        recal_table = "recal/{sample}.grp",
    log:
        "logs/recal/baserecalibrator_{sample}.log",
    params:
        extra = config["parameters"]["gatk"]["baserecalibrator"],
    resources:
        mem_mb=1024,
    wrapper:
        get_wrapper('gatk', 'baserecalibrator')


rule applybqsr:
    input:
        bam = "dedup/{sample}.rmdup.bam",
        ref = REFERENCE,
        dict = "resources/genome.dict",
        recal_table = rules.baserecalibrator.output.recal_table,
    output:
        bam = protected("recal/{sample}.bqsr.bam"),
    log:
        "logs/recal/applybqsr_{sample}.log",
    params:
        extra = config["parameters"]["gatk"]["applybqsr"],
    resources:
        mem_mb = 1024,
    wrapper:
        get_wrapper('gatk', 'applybqsr')


rule samtools_index:
    input:
        rules.applybqsr.output.bam,
    output:
        rules.applybqsr.output.bam + ".bai",
    log:
        "logs/recal/applybqsr_index_{sample}.log",
    wrapper:
        get_wrapper("samtools", "index")
