# rule fastp:
#     input:
#         reads = unpack(get_fastq),
#     output:
#         trimmed = ["trimmed/{sample}.clean.R1.fq.gz", "trimmed/{sample}.clean.R2.fq.gz"],
#         # Unpaired reads separately
#         # unpaired1 = "fastp/{sample}.u1.fastq",
#         # unpaired2 = "fastp/{sample}.u2.fastq",
#         # or in a single file
#         unpaired = "trimmed/{sample}.singletons.fastq",
#         # merged = "fastp/{sample}.merged.fastq",
#         failed = "trimmed/{sample}.failed.fastq",
#         html = "trimmed/report/{sample}.html",
#         json = "trimmed/report/{sample}.json"
#     log:
#         "logs/fastp_{sample}.log"
#     params:
#         extra = config['parameters']['fastp']['extra'],  # optional parameters
#         adapters = get_adapter('fastp')
#     threads: 2
#     wrapper:
#         get_wrapper('fastp')


rule trimmomatic:
    input:
        unpack(get_fastq),
    output:
        fq1 = "trimmed/{sample}.clean.R1.fq.gz",
        fq2 = "trimmed/{sample}.clean.R2.fq.gz",
        # reads where trimming entirely removed the mate
        fq1_unpaired = "trimmed/{sample}.R1.unpaired.fastq.gz",
        fq2_unpaired = "trimmed/{sample}.R2.unpaired.fastq.gz"
    log:
        "logs/trimmomatic_{sample}.log"
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

# rule cutadapt:
#     input:
#         ["reads/{sample}.1.fastq", "reads/{sample}.2.fastq"],
#     output:
#         fastq1 = "trimmed/{sample}.clean.R1.fq.gz",
#         fastq2 = "trimmed/{sample}.clean.R2.fq.gz",
#         qc = "trimmed/report/{sample}.qc.txt",
#     params:
#         # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
#         adapters = get_adapter('cutadapt'),
#         # https://cutadapt.readthedocs.io/en/stable/guide.html#
#         extra = "--minimum-length 1 -q 20",
#     log:
#         "logs/cutadapt/{sample}.log",
#     threads: 4  # set desired number of threads here
#     wrapper:
#         get_wrapper("cutadapt")


# rule fastqc:
#     input:
#         ["trimmed/{sample}.clean.R1.fq.gz", "trimmed/{sample}.clean.R2.fq.gz"]
#     output:
#         html = "qc/fastqc/{sample}.html",
#         zip = "qc/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
#     params: "--quiet"
#     log:
#         "logs/fastqc/{sample}.log"
#     threads: 1
#     wrapper:
#         get_wrapper("fastqc")