# 4. RealignerTargetCreator
rule realignertargetcreator:
    input:
        bam = "dedup/{sample}.rmdup.bam",
        bai = "dedup/{sample}.rmdup.bai",
        ref = REFERENCE,
        dict = "resources/genome.dict",
        known = KNOWN_SNP,
    output:
        intervals = "recal/{sample}.intervals",
    log:
        "logs/reacal/realignertargetcreator_{sample}.log",
    params:
        extra = "--defaultBaseQualities 20 --filter_reads_with_N_cigar",  # optional
    resources:
        mem_mb = 1024,
    threads: 16
    wrapper:
        get_wrapper('gatk3', 'realignertargetcreator')
# 
rule indelrealigner:
    input:
        bam = "dedup/{sample}.rmdup.bam",
        bai = "dedup/{sample}.rmdup.bai",
        ref = REFERENCE,
        dict = "resources/genome.dict",
        known = KNOWN_INDEL,
        target_intervals = "recal/{sample}.intervals",
    output:
        bam = "recal/{sample}.realigned.bam",
        bai = "recal/{sample}.realigned.bai",
    log:
        "logs/recal/indelrealigner_{sample}.log",
    params:
        extra="--defaultBaseQualities 20 --filter_reads_with_N_cigar",  # optional
    threads: 16
    resources:
        mem_mb=1024,
    wrapper:
        get_wrapper('gatk3', 'indelrealigner')
# BaseRecalibrator
rule baserecalibrator:
    input:
        bam = "dedup/{sample}.rmdup.bam",,
        ref = REFERENCE,
        dict = "resources/genome.dict",
        known = KNOWN,
    output:
        recal_table = "recal/{sample}.recal_data_table",
    log:
        "logs/recal/baserecalibrator_{sample}.log",
    params:
        extra="--defaultBaseQualities 20 --filter_reads_with_N_cigar",  # optional
    resources:
        mem_mb=1024,
    threads: 16
    wrapper:
        get_wrapper('gatk3', 'baserecalibrator')
# 
rule printreads:
    input:
        bam = "dedup/{sample}.rmdup.bam",
        bqsr = rules.baserecalibrator.output.recal_table,
        ref = REFERENCE,
        dict = "resources/genome.dict",
    output:
        bam = "recal/{sample}.bqsr.bam",
        bai = "recal/{sample}.bqsr.bai",
    log:
        "logs/recal/printreads_{sample}.log",
    params:
        extra = "--defaultBaseQualities 20 --filter_reads_with_N_cigar",  # optional
    resources:
        mem_mb = 1024,
    threads: 16
    wrapper:
        get_wrapper('gatk3', 'printreads')
