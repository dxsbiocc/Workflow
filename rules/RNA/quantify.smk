if QUANTIFY_TOOL == "rsem":
    if not QUANTIFY_INDEX:
        RSEM_INDEX = opj(OUTDIR, "quantity/index/reference")
        rule rsem_index:
            input:
                reference_genome = REFERENCE,
                gtf = GTF,
            output:
                seq = f"{RSEM_INDEX}.seq",
                grp = f"{RSEM_INDEX}.grp",
                ti = f"{RSEM_INDEX}.ti",
            params:
                extra = "--bowtie2",
            log:
                opj(OUTDIR, "logs/rsem/prepare-reference.log")
            threads: 8
            wrapper:
                get_wrapper("rsem", "prepare-reference")
    else:
        RSEM_INDEX = QUANTIFY_INDEX
    rule rsem_quant:
        input:
            bam = opj(OUTDIR, "dedup/{sample}/{sample}.rmdup.bam"),
            reference = multiext(RSEM_INDEX, ".grp", ".ti", ".seq"),
        output:
            results = directory(opj(OUTDIR, "quantity/{sample}/")),
            # genes_results = opj(OUTDIR, "quantity/{sample}/{sample}.genes.results"),
            # isoforms_results = opj(OUTDIR, "quantity/{sample}/{sample}.isoforms.results"),
        params:
            paired_end = PAIRED,
            mapping = MAPPING,
            extra = config["parameters"]["rsem"]["extra"],
        log:
            opj(OUTDIR, "logs/rsem/{sample}.log")
        threads: 10
        wrapper:
            get_wrapper("rsem", "calculate-expression")
            
elif QUANTIFY_TOOL == "salmon":
    if not QUANTIFY_INDEX:
        SALMON_INDEX = opj(OUTDIR, "quantity/index/reference")
        rule salmon_index:
            input:
                sequences = REFERENCE,
            output:
                directory(SALMON_INDEX),
            params:
                extra = "",
            log:
                opj(OUTDIR, "logs/salmon/index.log")
            threads: 8
            wrapper:
                get_wrapper("salmon", "index")
    else:
        SALMON_INDEX = QUANTIFY_INDEX
    rule salmon:
        input:
            # If you have multiple fastq files for a single sample (e.g. technical replicates)
            # use a list for r1 and r2.
            r1 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R1.fq.gz"),
            r2 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R2.fq.gz"),
            index = SALMON_INDEX,
        output:
            quant = directory(opj(OUTDIR, "quantity/{sample}")),
            # quant = opj(OUTDIR, "quantity/{sample}/{sample}.quant.sf"),
            lib = opj(OUTDIR, "quantity/{sample}/{sample}.lib_format_counts.json"),
        log:
            opj(OUTDIR, "logs/salmon/{sample}.log"),
        threads: 10
        params:
            # optional parameters
            libtype = config["parameters"]["salmon"]["libtype"],
            extra = config["parameters"]["salmon"]["extra"],
        wrapper:
            get_wrapper("salmon", "quant")
elif QUANTIFY_TOOL == "kallisto":
    if not QUANTIFY_INDEX:
        KALLISTO_INDEX = opj(OUTDIR, "quantity/index/reference")
        rule kallisto_index:
            input:
                fasta = REFERENCE,
            output:
                index = KALLISTO_INDEX,
            params:
                extra = "",
            log:
                opj(OUTDIR, "logs/kallisto/index.log")
            threads: 8
            wrapper:
                get_wrapper("kallisto", "index")
    else:
        KALLISTO_INDEX = QUANTIFY_INDEX
    rule kallisto:
        input:
            fastq = expand(opj(OUTDIR, "trimmed/{{sample}}/{{sample}}.clean.{run}.fq.gz"), run=RUN),
            index = KALLISTO_INDEX,
        output:
            directory(opj(OUTDIR, "quantity/{sample}/")),
        log:
            opj(OUTDIR, "logs/kallisto/{sample}.log"),
        threads: 10
        params:
            # optional parameters
            extra = config["parameters"]["kallisto"]["extra"],
        wrapper:
            get_wrapper("kallisto", "quant")
elif QUANTIFY_TOOL == "featurecounts":
    rule feature_counts:
        input:
            # list of sam or bam files
            samples = opj(OUTDIR, "dedup/{sample}/{sample}.rmdup.bam"),
            annotation = GTF,
            # optional input
            #chr_names="",           # implicitly sets the -A flag
            #fasta="genome.fasta"    # implicitly sets the -G flag
        output:
            directory(opj(OUTDIR, "quantity/{sample}"))
        threads: 10
        params:
            strand = config["parameters"]["featurecounts"]["strand"],  
            r_path = config["parameters"]["featurecounts"]["r_path"],
            extra = config["parameters"]["featurecounts"]["extra"],
            paired = PAIRED
        log:
            opj(OUTDIR, "logs/featurecounts/{sample}.log"),
        wrapper:
            get_wrapper("subread", "featurecounts")
else:
    raise ValueError("QUANTIFY_TOOL must be one of 'rsem', 'salmon', 'featurecounts', or 'kallisto'")

rule normalization:
    pass