rule hicfindrestsite:
    input:
        opj(OUTDIR, "genome/genome.{hap}.fasta")
    output:
        opj(OUTDIR, "scaffolding/restsite/{hap}.DpnII.bed")
    log:
        opj(OUTDIR, "logs/hicexplorer/hicfindrestsite_{hap}.log")
    params:
        pattern = ["GATC"],  # HindIII: 'AAGCTT', DpnII/MboI/Sau3AI: 'GATC', Arima: 'GATC', 'GANTC'
        extra = ""
    threads: 20
    wrapper:
        get_wrapper("hicexplorer", "pre-processing", "hicfindrestsite")

rule chromsize:
    input:
        opj(OUTDIR, "genome/genome.{hap}.fasta.fai")
    output:
        opj(OUTDIR, "genome/genome.{hap}.chrom.size")
    log:
        opj(OUTDIR, "logs/genome/{hap}.chromsize.log")
    shell:
        """
        awk -F'\t' '{{print $1, $2}}' {input} > {output}
        """

rule bwa_index:
    input:
        opj(OUTDIR, "genome/genome.{hap}.fasta"),
    output:
        idx = multiext(opj(OUTDIR, "genome/genome.{hap}.fasta"), ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        opj(OUTDIR, "logs/bwa_index/{hap}.log"),
    params:
        algorithm = "bwtsw",
    wrapper:
        get_wrapper("bwa", "index")

rule juicer:
    input:
        genome = opj(OUTDIR, "genome/genome.{hap}.fasta"),
        restriction_sites = rules.hicfindrestsite.output[0],
        chromsizes = rules.chromsize.output[0],
        fastq = expand(opj(OUTDIR, "trimmed/{hic}/{hic}.clean.{run}.fq.gz"), hic=DATA_DICT['HIC'], run=RUN),
        juicer = "/disk2/dxs/software/juicer"
    output:
        protected(directory(opj(OUTDIR, "scaffolding/juicer/{hap}")))
    log:
        opj(OUTDIR, "logs/juicer/build_{hap}.log")
    params:
        gname = "PlaTyr",
        restriction_type = 'DpnII',
    threads: 100
    wrapper:
        get_wrapper("juicer", "build")

rule threed_dna:
    input:
        fasta = opj(OUTDIR, "genome/genome.{hap}.fasta"),
        mnd = opj(OUTDIR, "scaffolding/juicer/{hap}/aligned/merged_nodups.txt"),
    output:
        opj(OUTDIR, "scaffolding/3d-dna/{hap}")
    log:
        opj(OUTDIR, "logs/3d-dna/{hap}.log")
    params:
        mode = 'haploid',    # haploid/diploid
        input_size = 15000,
        rounds = 2,          # Specifies number of iterative rounds for misjoin correction (default is 2)
        stage = '',          # Fast forward to later assembly steps, can be polish, split, seal, merge and finalize
    wrapper:
        get_wrapper("3d-dna")

rule bedtools_bamtobed:
    input:
        opj(OUTDIR, "scaffolding/juicer/{hap}/aligned/merged_dedup.bam"),
    output:
        opj(OUTDIR, "scaffolding/salsa2/{hap}.bed")
    log:
        opj(OUTDIR, "logs/salsa2/bamtobed_{hap}.log")
    params:
        extra = "-bedpe",  # optional
    wrapper:
        get_wrapper("bedtools", "bamtobed")

rule salsa2:
    input:
        fas = opj(OUTDIR, "genome/genome.{hap}.fasta"),
        fai = opj(OUTDIR, "genome/genome.{hap}.fasta.fai"),
        bed = opj(OUTDIR, "scaffolding/salsa2/{hap}.bed"),
    output:
        directory(opj(OUTDIR, "scaffolding/salsa2/{hap}"))
    log:
        opj(OUTDIR, "logs/salsa2/{hap}.log"),
    params:
        enzyme = "GATC",  # optional
        extra = "--clean yes -m yes",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        get_wrapper("salsa2")
