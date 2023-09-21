rule gfatools_stat:
    input:
        "assembly/genome.{suffix}.gfa",
    output:
        "genome/genome.{suffix}.stat",
    log:
        "logs/gfatools/{suffix}.stat.log",
    params:
        command = "stat",
    wrapper:
        get_wrapper("gfatools")


rule gfatools_gfa2fa:
    input:
        "assembly/genome.{suffix}.gfa",
    output:
        "genome/genome.{suffix}.fasta",
    log:
        "logs/gfatools/{suffix}.gfa2fa.log",
    params:
        command = "gfa2fa",
	    extra = "",
    wrapper:
        get_wrapper("gfatools")


rule samtools_index:
    input:
        "genome/genome.{suffix}.fasta",
    output:
        "genome/genome.{suffix}.fasta.fai",
    log:
        "logs/index/{suffix}.log",
    params:
        extra = "",  # optional params string
    wrapper:
        get_wrapper("samtools", "faidx")

rule seqkit_stats:
    input:
        fastx = "genome/genome.{suffix}.fasta",
    output:
        stats="qc/stats/{suffix}.tsv",
    log:
        "logs/seqkit_stats/{suffix}.log",
    params:
        command = "stats",
        extra = "--all --tabular",
    threads: 20
    wrapper:
        get_wrapper("seqkit")


rule meryl_count:
    input:
        fasta = expand("trimmed/{hifi}/{hifi}.filt.fastq.gz", hifi=DATA_DICT['HIFI']),
    output:
        directory("merqury/genome/"),
    log:
        "logs/merqury/meryl_count.log",
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
        hifi = expand("genome/genome.{suffix}.fasta", suffix=["hap1.p_ctg", "hap2.p_ctg"]),
        db = rules.meryl_count.output,
    output:
        directory("merqury/diploid")
    log:
        "logs/merqury/diploid.log",
    params:
        output_prefix = "merqury",
        extra = "",
    threads: 100
    wrapper:
        get_wrapper("merqury")


rule busco_hap:
    input:
        fasta = "genome/genome.{suffix}.fasta",
    output:
        directory("busco/{suffix}"),
    log:
        "logs/qc/busco_{suffix}.log",
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