if MAPPING == "bwa":
    # DNA-seq
    rule bwa:
        input:
            reads = expand(opj(OUTDIR, "trimmed/{{sample}}/{{sample}}.clean.{run}.fq.gz"), run=RUN),
            index = multiext(
                INDEX, 
                ".amb", ".ann", ".bwt", ".pac", ".sa"
            )
        output:
            opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam"),
        log:
            opj(OUTDIR, "logs/mapping/bwa_{sample}.log"),
        params:
            extra = config['parameters']['bwa']['extra'],
            sort = config['parameters']['bwa']['sort'],  # Can be 'none', 'samtools' or 'picard'.
            sort_order = config['parameters']['bwa']['sort_order'],  # Can be 'queryname' or 'coordinate'.
            sort_extra = config['parameters']['bwa']['sort_extra'],  # Extra args for samtools/picard.
        threads: config['parameters']['bwa']['threads']
        wrapper:
            get_wrapper("bwa", "mem")
elif MAPPING == "bwa-mem2":
    rule bwa_mem2_mem:
        input:
            reads = expand(opj(OUTDIR, "trimmed/{{sample}}/{{sample}}.clean.{run}.fq.gz"), run=RUN),
            # Index can be a list of (all) files created by bwa, or one of them
            idx = multiext(INDEX, ".amb", ".ann", ".bwt.2bit.64", ".pac"),
        output:
            opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam"),
        log:
            opj(OUTDIR, "logs/mapping/bwa2_{sample}.log"),
        params:
            extra = config['parameters']['bwa']['extra'],
            sort = config['parameters']['bwa']['sort'],  # Can be 'none', 'samtools' or 'picard'.
            sort_order = config['parameters']['bwa']['sort_order'],  # Can be 'queryname' or 'coordinate'.
            sort_extra = config['parameters']['bwa']['sort_extra'],  # Extra args for samtools/picard.
        threads: config['parameters']['bwa']['threads']
        wrapper:
            get_wrapper("bwa-mem2", "mem")
elif MAPPING == "bowtie2":
    # RNA-seq, short-read
    rule bowtie2:
        input:
            reads = expand(opj(OUTDIR, "trimmed/{{sample}}/{{sample}}.clean.{run}.fq.gz"), run=RUN),
            index = multiext(
                INDEX,
                ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
            )
        output:
            opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam"),
        log:
            opj(OUTDIR, "logs/mapped/bowtie2_{sample}.summary"),
        params:
            extra = config['parameters']['bowtie2']['extra'],  # optional parameters
            sort = config['parameters']['bowtie2']['sort']
        threads: config['parameters']['bowtie2']['threads'] # Use at least two threads
        wrapper:
            get_wrapper("bowtie2", "align")
elif MAPPING == "hisat2":
    rule hisat2:
        input:
            reads = expand(opj(OUTDIR, "trimmed/{{sample}}/{{sample}}.clean.{run}.fq.gz"), run=RUN),
            index = multiext(
                INDEX,
                ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2", ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"
            )
        output:
            opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam")
        log:
            opj(OUTDIR, "logs/mapped/hisat2_{sample}.log")
        params:
            extra = config['parameters']['hisat2']['extra'],
            sort = config['parameters']['hisat2']['sort'],              # Can be 'none', 'samtools' or 'picard'.
            sort_order = config['parameters']['hisat2']['sort_order'],  # Can be 'queryname' or 'coordinate'.
            sort_extra = config['parameters']['hisat2']['sort_extra'],
        threads: config['parameters']['hisat2']['threads']
        wrapper:
            get_wrapper("hisat2", "align")
elif MAPPING == "star":
    rule star:
        input:
            reads = expand(opj(OUTDIR, "trimmed/{{sample}}/{{sample}}.clean.{run}.fq.gz"), run=RUN),
            # path to STAR reference genome index
            index = INDEX,
        output:
            # see STAR manual for additional output files
            aln = opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam"),
            unmapped = expand(opj(OUTDIR, "mapped/{{sample}}/{{sample}}.unmapped.{run}.fastq"), run=RUN),
        log:
            opj(OUTDIR, "logs/mapped/star_{sample}.log"),
        params:
            # optional parameters
            extra = config['parameters']['star']['extra'] + \
                    " --quantTranscriptomeBan {}".format("IndelSoftclipSingleend") if QUANTIFY_TOOL == "rsem" else "" + \
                    " --outSAMattrRGline ID:{sample} SM:{sample}" + \
                    " --sjdbGTFfile {}".format(GTF) + \
                    " --sjdbOverhang {}".format(READ_LENGTH - 1),
            fusion = config['control'].get('fusion', '').lower() if config['control'].get('merge') else ""
        threads: config['parameters']['star']['threads']
        wrapper:
            get_wrapper("star", "align")
elif MAPPING == "minimap2":
    # long read mapping
    rule minimap2:
        input:
            query = expand(opj(OUTDIR, "trimmed/{{sample}}/{{sample}}.clean.{run}.fq.gz"), run=RUN),
            target = INDEX,  # can be either genome index or genome fasta
        output:
            opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam"),
        log:
            opj(OUTDIR, "logs/mapped/minimap2_{sample}.log"),
        params:
            extra = config['parameters']['minimap2']['extra'],            # optional
            sorting = config['parameters']['minimap2']['sorting'],        # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
            sort_extra = config['parameters']['minimap2']['sort_extra'],  # optional: extra arguments for samtools/picard
        threads: config['parameters']['minimap2']['threads']
        wrapper:
            get_wrapper("minimap2", "aligner")
else:
    raise ValueError(f'This alignment tool: {MAPPING} not support!')