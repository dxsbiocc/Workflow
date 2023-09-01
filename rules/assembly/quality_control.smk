rule gfastats:
    input:
        "assembly/genome.{suffix}.p_ctg.gfa",
    output:
        seq = "genome/genome.{suffix}.fasta",        # output format (fasta/fastq/gfa)
        summary = "genome/genome.{suffix}.assembly_summary",
    log:
        "logs/genome/{suffix}_gfastats.log",
    params:
        extra = "",
        agpfile = "",
        include_bed = "",
        exclude_bed = "",
        instructions = "",
    threads: 10
    wrapper:
        get_wrapper("gfastats")

rule run_busco:
    input:
        rules.gfastats.output.seq,
    output:
        out_dir = directory("busco/{suffix}"),
        dataset_dir = directory("busco/downloads"),
    log:
        "logs/qc/busco_{suffix}.log",
    params:
        mode = "genome",
        # optional parameters
        extra = "--auto-lineage-euk",
    threads: 10
    wrapper:
        get_wrapper("busco")

rule DNA_realign:  # ratio > 95%
    pass

rule RNA_realign:  # ratio > 80%
    pass

rule GC_depth:     # single peak
    pass

rule remove_contamination:
    pass

rule SNP_detection:
    pass