if QUANTIFY_TOOL == "rsem":
    if not RSEM_INDEX:
        RSEM_INDEX = opj(OUTDIR, "quantity/index/reference")
        rule rsem_index:
            input:
                reference_genome = REFERENCE,
            output:
                seq = f"{RSEM_INDEX}.seq",
                grp = f"{RSEM_INDEX}.grp",
                ti = f"{RSEM_INDEX}.ti",
            params:
                extra = "",
            log:
                opj(OUTDIR, "logs/rsem/prepare-reference.log")
            threads: 8
            wrapper:
                get_wrapper("rsem", "prepare-reference")
    rule rsem_quant:
        input:
            bam = opj(OUTDIR, "dedup/{sample}/{sample}.rmdup.bam"),
            reference = multiext(RSEM_INDEX, ".grp", ".ti", ".transcripts.fa", ".seq", ".idx.fa", ".n2g.idx.fa"),
        output:
            genes_results = opj(OUTDIR, "quantity/{sample}/{sample}.genes.results"),
            isoforms_results = opj(OUTDIR, "quantity/{sample}/{sample}.isoforms.results"),
        params:
            paired_end = PAIRED,
            extra = config["parameters"]["rsem"]["extra"],
        log:
            opj(OUTDIR, "logs/rsem/{sample}.log")
        wrapper:
            get_wrapper("rsem", "calculate-expression")
            
elif QUANTIFY_TOOL == "salmon":
    if not SALMON_INDEX:
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
    rule salmon:
        input:
            # If you have multiple fastq files for a single sample (e.g. technical replicates)
            # use a list for r1 and r2.
            r1 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R1.fq.gz"),
            r2 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R2.fq.gz"),
            index = SALMON_INDEX,
        output:
            quant = opj(OUTDIR, "quantity/{sample}/{sample}.quant.sf"),
            lib = opj(OUTDIR, "quantity/{sample}/{sample}.lib_format_counts.json"),
        log:
            opj(OUTDIR, "logs/salmon/{sample}.log"),
        params:
            # optional parameters
            libtype = config["parameters"]["salmon"]["libtype"],
            extra = config["parameters"]["salmon"]["extra"],
        wrapper:
            get_wrapper("salmon", "quant")
elif QUANTIFY_TOOL == "kallisto":
    if not KALLISTO_INDEX:
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
    rule kallisto:
        input:
            fastq = expand(opj(OUTDIR, "trimmed/{{sample}}/{{sample}}.clean.{run}.fq.gz"), run=RUN),
            index = KALLISTO_INDEX,
        output:
            directory(opj(OUTDIR, "quantity/{sample}/")),
        log:
            opj(OUTDIR, "logs/kallisto/{sample}.log"),
        params:
            # optional parameters
            extra = config["parameters"]["kallisto"]["extra"],
        wrapper:
            get_wrapper("kallisto", "quant")
else:
    raise ValueError("QUANTIFY_TOOL must be one of 'rsem', 'salmon', or 'kallisto'")

rule normalization:
    pass