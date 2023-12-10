# Run minimap2 to align pacbio data and generate paf files, 
# then calculate read depth histogram and base-level read depth
rule minimap2_paf:
    input:
        target = opj(OUTDIR, "genome/genome.p_ctg.fasta"),  # can be either genome index or genome fasta
        query = opj(OUTDIR, "trimmed/{hifi}/{hifi}.filt.fastq.gz"),
    output:
        opj(OUTDIR, "dedups/{hifi}.paf.gz"),
    log:
        opj(OUTDIR, "logs/minimap2/{hifi}.log"),
    params:
        extra = "-x map-pb",  # optional
        sorting = "coordinate",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra = "",  # optional: extra arguments for samtools/picard
    threads: 200
    wrapper:
        get_wrapper("minimap2", "aligner")

rule purge_dups_pbcstat:
    input:
        expand(opj(OUTDIR, "dedups/{hifi}.paf.gz"), hifi=DATA_DICT["hifi"]),
    output:
        cov = opj(OUTDIR, "dedups/purge_dups/pbcstat/genome.cov"),
        stat = opj(OUTDIR, "dedups/purge_dups/pbcstat/genome.stat"),
    log:
        opj(OUTDIR, "logs/purge_dups/pbcstat.log"),
    params:
        extra = "",
    threads: 30
    wrapper:
        get_wrapper("purge_dups", "pbcstat")

rule purge_dups_calcuts:
    input:
        rules.purge_dups_pbcstat.output.stat,
    output:
        opj(OUTDIR, "dedups/purge_dups/calcuts/genome.cutoffs"),
    log:
        opj(OUTDIR, "logs/purge_dups/calcuts.log"),
    params:
        extra = "-l 2 -m 4 -u 8",
    threads: 30
    wrapper:
        get_wrapper("purge_dups", "calcuts")

# Split an assembly and do a self-self alignment
rule purge_dups_split_fa:
    input:
        opj(OUTDIR, "genome/genome.p_ctg.fasta"),
    output:
        opj(OUTDIR, "dedups/purge_dups/split_fa/primary.split"),
    log:
        opj(OUTDIR, "logs/purge_dups/split_fa.log"),
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
        opj(OUTDIR, "dedups/primary.split.self.paf.gz"),
    log:
        opj(OUTDIR, "logs/minimap2/split.self.log"),
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
        opj(OUTDIR, "dedups/purge_dups/purge_dups/genome.bed"),
    log:
        opj(OUTDIR, "logs/purge_dups/purge_dups.log"),
    params:
        extra = "-2",
    threads: 50
    wrapper:
        get_wrapper("purge_dups", "purge_dups")

# Get purged primary and haplotig sequences from draft assembly
rule purge_dups_get_seqs:
    input:
        fas = opj(OUTDIR, "genome/genome.p_ctg.fasta"),
        bed = rules.purge_dups.output,
    output:
        hap = opj(OUTDIR, "dedups/get_seqs.hap.fasta"),
        purged = opj(OUTDIR, "dedups/get_seqs.purged.fasta"),
    log:
        opj(OUTDIR, "logs/purge_dups/get_seqs.log"),
    params:
        extra = "",
    threads: 20
    wrapper:
        get_wrapper("purge_dups", "get_seqs")