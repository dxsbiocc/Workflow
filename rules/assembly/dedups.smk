# Run minimap2 to align pacbio data and generate paf files, 
# then calculate read depth histogram and base-level read depth
rule minimap2_paf:
    input:
        target = "genome/genome.p_ctg.fasta",  # can be either genome index or genome fasta
        query = "trimmed/{hifi}/{hifi}.filt.fastq.gz",
    output:
        "dedups/{hifi}.paf.gz",
    log:
        "logs/minimap2/{hifi}.log",
    params:
        extra = "-x map-pb",  # optional
        sorting = "coordinate",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra = "",  # optional: extra arguments for samtools/picard
    threads: 200
    wrapper:
        get_wrapper("minimap2", "aligner")

rule purge_dups_pbcstat:
    input:
        expand("dedups/{hifi}.paf.gz", hifi=DATA_DICT["hifi"]),
    output:
        cov = "dedups/purge_dups/pbcstat/genome.cov",
        stat = "dedups/purge_dups/pbcstat/genome.stat",
    log:
        "logs/purge_dups/pbcstat.log",
    params:
        extra = "",
    threads: 30
    wrapper:
        get_wrapper("purge_dups", "pbcstat")

rule purge_dups_calcuts:
    input:
        rules.purge_dups_pbcstat.output.stat,
    output:
        "dedups/purge_dups/calcuts/genome.cutoffs",
    log:
        "logs/purge_dups/calcuts.log",
    params:
        extra = "-l 2 -m 4 -u 8",
    threads: 30
    wrapper:
        get_wrapper("purge_dups", "calcuts")

# Split an assembly and do a self-self alignment
rule purge_dups_split_fa:
    input:
        "genome/genome.p_ctg.fasta",
    output:
        "dedups/purge_dups/split_fa/primary.split",
    log:
        "logs/purge_dups/split_fa.log",
    params:
        extra = "",
    threads: 30
    wrapper:
        get_wrapper("purge_dups", "split_fa")

rule minimap2_self:
    input:
        target = rules.purge_dups_split_fa.output,  # can be either genome index or genome fasta
        query = rules.purge_dups_split_fa.output,
    output:
        "dedups/primary.split.self.paf.gz",
    log:
        "logs/minimap2/split.self.log",
    params:
        extra = "-xasm5 -DP",  # optional
        sorting = "coordinate",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra = "",  # optional: extra arguments for samtools/picard
    threads: 200
    wrapper:
        get_wrapper("minimap2", "aligner")

# Purge haplotigs and overlaps
rule purge_dups:
    input:
        paf = rules.minimap2_self.output,
        cov = rules.purge_dups_pbcstat.output.cov,
        cutoff = rules.purge_dups_calcuts.output,
    output:
        "dedups/purge_dups/purge_dups/genome.bed",
    log:
        "logs/purge_dups/purge_dups.log",
    params:
        extra = "-2",
    threads: 50
    wrapper:
        get_wrapper("purge_dups", "purge_dups")

# Get purged primary and haplotig sequences from draft assembly
rule purge_dups_get_seqs:
    input:
        fas = "genome/genome.p_ctg.fasta",
        bed = rules.purge_dups.output,
    output:
        hap = "dedups/get_seqs.hap.fasta",
        purged = "dedups/get_seqs.purged.fasta",
    log:
        "logs/purge_dups/get_seqs.log",
    params:
        extra = "",
    threads: 20
    wrapper:
        get_wrapper("purge_dups", "get_seqs")