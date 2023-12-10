rule gfatools_stat:
    input:
        opj(OUTDIR, "assembly/genome.{suffix}.gfa"),
    output:
        opj(OUTDIR, "genome/genome.{suffix}.stat"),
    log:
        opj(OUTDIR, "logs/gfatools/{suffix}.stat.log"),
    params:
        command = "stat",
    wrapper:
        get_wrapper("gfatools")


rule gfatools_gfa2fa:
    input:
        opj(OUTDIR, "assembly/genome.{suffix}.gfa"),
    output:
        opj(OUTDIR, "genome/genome.{suffix}.fasta"),
    log:
        opj(OUTDIR, "logs/gfatools/{suffix}.gfa2fa.log"),
    params:
        command = "gfa2fa",
	    extra = "",
    wrapper:
        get_wrapper("gfatools")


rule samtools_index:
    input:
        opj(OUTDIR, "genome/genome.{suffix}.fasta"),
    output:
        opj(OUTDIR, "genome/genome.{suffix}.fasta.fai"),
    log:
        opj(OUTDIR, "logs/index/{suffix}.log"),
    params:
        extra = "",  # optional params string
    wrapper:
        get_wrapper("samtools", "faidx")

rule seqkit_stats:
    input:
        fastx = opj(OUTDIR, "genome/genome.{suffix}.fasta"),
    output:
        stats = opj(OUTDIR, "qc/stats/{suffix}.tsv"),
    log:
        opj(OUTDIR, "logs/seqkit_stats/{suffix}.log"),
    params:
        command = "stats",
        extra = "--all --tabular",
    threads: 20
    wrapper:
        get_wrapper("seqkit")


rule meryl_count:
    input:
        fasta = expand(opj(OUTDIR, "trimmed/{hifi}/{hifi}.filt.fastq.gz"), hifi=DATA_DICT['HIFI']),
    output:
        directory(opj(OUTDIR, "merqury/genome/")),
    log:
        opj(OUTDIR, "logs/merqury/meryl_count.log"),
    params:
        command = "count",
        extra = "k=31",
    threads: 200
    resources:
        mem_mb = 504800,
    wrapper:
        get_wrapper("meryl", "count")


rule merqury_diploid:
    input:
        hifi = expand(opj(OUTDIR, "genome/genome.{suffix}.fasta"), suffix=["hap1.p_ctg", "hap2.p_ctg"]),
        db = rules.meryl_count.output,
    output:
        directory(opj(OUTDIR, "merqury/diploid"))
    log:
        opj(OUTDIR, "logs/merqury/diploid.log"),
    params:
        output_prefix = "merqury",
        extra = "",
    threads: 100
    wrapper:
        get_wrapper("merqury")


rule busco_hap:
    input:
        fasta = opj(OUTDIR, "genome/genome.{suffix}.fasta"),
    output:
        directory(opj(OUTDIR, "busco/{suffix}")),
    log:
        opj(OUTDIR, "logs/qc/busco_{suffix}.log"),
    params:
        dataset_dir = "vertebrata_odb10"
        mode = "genome",
        lineage = "vertebrata_odb10",
        # optional parameters
        # extra = "--auto-lineage-euk",
    threads: 10
    wrapper:
        get_wrapper("busco")

# rule DNA_realign:  # ratio > 95%
#     pass

# rule RNA_realign:  # ratio > 80%
#     pass

# rule GC_depth:     # single peak
#     pass

# rule remove_contamination:
#     pass

# rule SNP_detection:
#     pass