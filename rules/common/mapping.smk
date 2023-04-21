# DNA-seq
rule bwa:
    input:
        reads = expand("trimmed/{{sample}}.clean.{run}.fq.gz", run=RUN),
        index = multiext(
            INDEX, 
            ".amb", ".ann", ".bwt", ".pac", ".sa"
        ),
    output:
        "mapped/{sample}.sorted.bam",
    log:
        "logs/mapping/bwa_{sample}.log",
    params:
        extra = config['parameters']['bwa']['extra'],
        sort = config['parameters']['bwa']['sort'],  # Can be 'none', 'samtools' or 'picard'.
        sort_order = config['parameters']['bwa']['sort_order'],  # Can be 'queryname' or 'coordinate'.
        sort_extra = config['parameters']['bwa']['sort_extra'],  # Extra args for samtools/picard.
    threads: 2
    wrapper:
        get_wrapper("bwa", "mem")

# RNA-seq
rule bowtie2:
    input:
        reads = expand("trimmed/{{sample}}.clean.{run}.fq.gz", run=RUN),
        index = multiext(
            INDEX,
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        ),
    output:
        "mapped/{sample}.sorted.bam",
    log:
        "logs/mapped/bowtie2_{sample}.summary",
    params:
        extra = config['parameters']['bowtie2']['extra'],  # optional parameters
        sort = config['parameters']['bowtie2']['sort']
    threads: 8  # Use at least two threads
    wrapper:
        get_wrapper("bowtie2", "align")

rule hisat2:
    input:
        reads = expand("trimmed/{{sample}}.clean.{run}.fq.gz", run=RUN),
        index = multiext(
            INDEX,
            ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2", ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"
            )
    output:
        "mapped/{sample}.sorted.bam"
    log:
        "logs/mapped/hisat2_{sample}.log"
    params:
        extra = config['parameters']['hisat2']['extra'],
    threads: 2
    wrapper:
        get_wrapper("hisat2", "align")

rule star:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fastq1 = "trimmed/{sample}.clean.R1.fq.gz",
        # paired end reads needs to be ordered so each item in the two lists match
        fastq2 = "trimmed/{sample}.clean.R2.fq.gz",  #optional
        # path to STAR reference genome index
        index = INDEX,
    output:
        # see STAR manual for additional output files
        aln = "mapped/{sample}.sorted.bam",
        log = "mapped/Log.out",
        sj = "mapped/SJ.out.tab",
    log:
        "logs/mapped/star_{sample}.log",
    params:
        # optional parameters
        extra = config['parameters']['star']['extra'] + " --outSAMattrRGline ID:{sample} SM:{sample}",
    threads: 8
    wrapper:
        get_wrapper("star", "align")

# long read mapping
rule minimap2:
    input:
        target = INDEX,  # can be either genome index or genome fasta
        query =expand("trimmed/{{sample}}.clean.{run}.fq.gz", run=RUN),
    output:
        "mapped/{sample}.sorted.bam",
    log:
        "logs/mapped/minimap2_{sample}.log",
    params:
        extra = config['parameters']['minimap2']['extra'],           # optional
        sorting = config['parameters']['minimap2']['sorting'],        # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra = config['parameters']['minimap2']['sort_extra'],  # optional: extra arguments for samtools/picard
    threads: 3
    wrapper:
        get_wrapper("minimap2", "aligner")