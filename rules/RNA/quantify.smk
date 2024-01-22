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
                extra = "",
                mapping = MAPPING,
            log:
                opj(OUTDIR, "logs/rsem/prepare-reference.log")
            threads: 8
            wrapper:
                get_wrapper("rsem", "prepare-reference")
    else:
        RSEM_INDEX = QUANTIFY_INDEX
    rule rsem_quant:
        input:
            bam = opj(OUTDIR, "mapped/{sample}/Aligned.toTranscriptome.out.bam"),
            fq_one = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R1.fq.gz"),
            fq_two = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R2.fq.gz") if PAIRED else "",
            reference = multiext(RSEM_INDEX, ".grp", ".ti", ".seq"),
        output:
            quant = opj(OUTDIR, "quantity/{sample}/{sample}.quant")
        params:
            paired_end = PAIRED,
            mapping = MAPPING, 
            extra = config["parameters"]["rsem"]["extra"],
            mode = config["parameters"]["rsem"]["mode"],
        log:
            opj(OUTDIR, "logs/rsem/quant_{sample}.log")
        threads: 10
        wrapper:
            get_wrapper("rsem", "calculate-expression")
            
elif QUANTIFY_TOOL == "salmon":
    if not QUANTIFY_INDEX:
        SALMON_INDEX = opj(OUTDIR, "quantity/index/reference")
        rule salmon_index:
            input:
                sequences = REFERENCE if not RNA else RNA,
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
            r1 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R1.fq.gz"),
            r2 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R2.fq.gz"),
            index = SALMON_INDEX,
            # aln = opj(OUTDIR, "mapped/{sample}/Aligned.toTranscriptome.out.bam"),
            # targeted = RNA,
        output:
            quant = opj(OUTDIR, "quantity/{sample}/{sample}.quant")
            # quant = opj(OUTDIR, "quantity/{sample}/{sample}.quant.sf"),
            # lib = opj(OUTDIR, "quantity/{sample}/{sample}.lib_format_counts.json"),
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
                fasta = REFERENCE if not RNA else RNA,
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
            # directory(opj(OUTDIR, "quantity/{sample}/")),
            quant = opj(OUTDIR, "quantity/{sample}/{sample}.quant")
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
            samples = opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam"),
            annotation = GTF,
        output:
            quant = opj(OUTDIR, "quantity/{sample}/{sample}.quant")
        threads: 10
        params:
            extra = config["parameters"]["featurecounts"]["extra"],
            strand = config["parameters"]["featurecounts"]["strand"],  
            r_path = config["parameters"]["featurecounts"]["r_path"],
            paired = PAIRED
        log:
            opj(OUTDIR, "logs/featurecounts/{sample}.log"),
        wrapper:
            get_wrapper("subread", "featurecounts")
elif QUANTIFY_TOOL == "htseq":
    rule htseq:
        input:
            bam = opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam"),
            anno = GTF
        output:
            quant = opj(OUTDIR, "quantity/{sample}/{sample}.quant")
        threads: 10
        params:
            extra = config["parameters"]["htseq"]["extra"],
            mode = "intersection-nonempty",  # union, intersection-strict and intersection-nonempty
            stranded = "no",  # yes/no/reverse
            order = "pos",    # name or pos
        log:
            opj(OUTDIR, "logs/htseq/{sample}.log")
        wrapper:
            get_wrapper("htseq")
else:
    raise ValueError("QUANTIFY_TOOL must be one of 'rsem', 'salmon', 'featurecounts', or 'kallisto'")

rule merge_expression:
    input:
        quant = expand(opj(OUTDIR, "quantity/{sample}/{sample}.quant"), sample=SAMPLES),
    output:
        exp = opj(OUTDIR, "expression/expression.tsv"),
    params:
        quantifier = QUANTIFY_TOOL,
    log:
        opj(OUTDIR, "logs/merge_counts.log")
    threads: 1
    wrapper:
        get_wrapper("scripts", "Python", "merge_expression")