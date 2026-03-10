rule metawrap_assembly:
    input:
        r1 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R1.fq.gz"),
        r2 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R2.fq.gz"),
    output:
        opj(OUTDIR, "metawrap/assembly/{sample}/final.contigs.fasta"),
    log:
        opj(OUTDIR, "logs/metawrap/assembly/{sample}.log"),
    params:
        assembly_method = config['parameters']['metawrap']['assembly_method'],
        extra = config['parameters']['metawrap']['assembly_extra'],
    threads: config['parameters']['metawrap']['assembly_threads'],
    resources:
        mem_mb = config['parameters']['metawrap']['assembly_memory'] * 1024,
    wrapper:
        get_wrapper("metawrap", "assembly")

rule metawrap_binning:
    input:
        assembly = opj(OUTDIR, config['parameters']['metawrap']['assembly_method'], "{sample}/final.contigs.fasta") \
            if BINNING_WITHOUT_ASSEMBLY else rules.metawrap_assembly.output,
        fastq = expand(opj(OUTDIR, "trimmed/{{sample}}/{{sample}}.clean.{run}.fq.gz"), run=RUN),
    output:
        directory(opj(OUTDIR, "metawrap/binning/{sample}")),
    log:
        opj(OUTDIR, "logs/metawrap/binning/{sample}.log"),
    params:
        extra = config['parameters']['metawrap']['binning_extra'],
        methods = config['parameters']['metawrap']['binning_method'],
    threads: config['parameters']['metawrap']['binning_threads'],
    resources:
        mem_mb = config['parameters']['metawrap']['binning_memory'] * 1024,
    wrapper:
        get_wrapper("metawrap", "binning")

rule metawrap_bins_refine:
    input:
        bins_dir = rules.metawrap_binning.output,
    output:
        directory(opj(OUTDIR, "metawrap/bins_refine/{sample}")),
    log:
        opj(OUTDIR, "logs/metawrap/bins_refine/{sample}.log"),
    params:
        extra = config['parameters']['metawrap']['refinement_extra'],
        bins = config['parameters']['metawrap']['binning_method'],
        checkm_db = config['parameters']['metawrap'].get('checkm_db', ""),
    threads: config['parameters']['metawrap']['refinement_threads'],
    resources:
        mem_mb = config['parameters']['metawrap']['refinement_memory'] * 1024,
    wrapper:
        get_wrapper("metawrap", "refinement")

rule drep_input:
    input:
        expand(opj(OUTDIR, "metawrap/bins_refine/{sample}"), sample=SAMPLES),
    output:
        directory(opj(OUTDIR, "drep/input/")),
    log:
        opj(OUTDIR, "logs/drep/input.log"),
    run:
        for sample in SAMPLES:
            shell(f"ln -s {opj(OUTDIR, "metawrap/bins_refine/{sample}")} {opj(OUTDIR, "drep/input/{sample}")}")

rule drep:
    input:
        genomes = rules.drep_input.output,
    output:
        directory(opj(OUTDIR, "drep/{sample}")),
    log:
        opj(OUTDIR, "logs/drep/{sample}.log"),
    params:
        extra = config['parameters']['drep']['extra'],
        command = "dereplicate",
    threads: config['parameters']['drep']['threads'],
    wrapper:
        get_wrapper("drep")

rule coverm_genome:
    input:
        reads = expand(opj(OUTDIR, "trimmed/{{sample}}/{{sample}}.clean.{run}.fq.gz"), run=RUN),
        fasta_dir = rules.drep.output,
    output:
        opj(OUTDIR, "coverm/{sample}.abundance.tsv"),
    log:
        opj(OUTDIR, "logs/coverm/{sample}.log"),
    params:
        extra = config['parameters']['coverm']['extra'],
        command = "genome",
    threads: config['parameters']['coverm']['threads'],
    wrapper:
        get_wrapper("coverm")

rule coverm_merge:
    input:
        expand(opj(OUTDIR, "coverm/{sample}.abundance.tsv"), sample=SAMPLES),
    output:
        opj(OUTDIR, "coverm/abundance.tsv"),
    log:
        opj(OUTDIR, "logs/coverm/abundance.log"),
    params:
        extra = "--file_name tsv",
    wrapper:
        get_wrapper("humann", "join")