rule fastp:
    input:
        reads = get_fastq,
    output:
        trimmed = expand("trimmed/{{sample}}.clean.{run}.fq.gz", run=RUN),
        # or in a single file
        unpaired = "trimmed/{sample}.singletons.fastq",
        # merged = "fastp/{sample}.merged.fastq",
        failed = "trimmed/{sample}.failed.fastq",
        html = "trimmed/report/{sample}.html",
        json = "trimmed/report/{sample}.json"
    log:
        "logs/trimmed/fastp_{sample}.log"
    params:
        extra = config['parameters']['fastp']['extra'],  # optional parameters
        adapters = get_adapter('fastp')
    threads: 2
    wrapper:
        get_wrapper('fastp')

rule trimmomatic:
    input:
        get_fastq,
    output:
        trimmed = expand("trimmed/{{sample}}.clean.{run}.fq.gz", run=RUN),
        # reads where trimming entirely removed the mate
        unpaired = expand("trimmed/{{sample}}.unpaired.{run}.fq.gz", run=RUN),
    log:
        "logs/trimmed/trimmomatic_{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer = config['parameters']['trimmomatic']['trimmer'],
        # get adapters file
        extra = get_adapter('trimmomatic'),
        compression_level = config['parameters']['trimmomatic']['compression_level']
    threads:
        32
    resources:
        mem_mb = config['parameters']['trimmomatic']['mem']
    wrapper:
        get_wrapper("trimmomatic")

rule cutadapt:
    input:
        get_fastq,
    output:
        trimmed = expand("trimmed/{{sample}}.clean.{run}.fq.gz", run=RUN),
        qc = "trimmed/report/{sample}.qc.txt",
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters = get_adapter('cutadapt'),
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        extra = config['parameters']['cutadapt']['extra'],
    log:
        "logs/trimmed/cutadapt_{sample}.log",
    threads: 4  # set desired number of threads here
    wrapper:
        get_wrapper("cutadapt")

rule trim_galore:
    input:
        get_fastq,
    output:
        trimmed = expand("trimmed/{{sample}}.clean.{run}.fq.gz", run=RUN),
        report = expand("trimmed/{{sample}}.trimming_report.{run}.fq.gz", run=RUN),
    params:
        extra = config['parameters']['trim_galore']['extra'],
    log:
        "logs/trimmed/trim_galore_{sample}.log",
    wrapper:
        get_wrapper("trim_galore")