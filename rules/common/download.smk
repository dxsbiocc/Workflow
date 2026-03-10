rule metadata:
    output:
        opj(OUTDIR, "metadata.tsv")
    log:
        opj(OUTDIR, "logs/sratools/metadata.log")
    wrapper:
        get_wrapper('pysradb', 'metadata')

rule prefecth:
    output:
        output_dir = opj(OUTDIR, "SRA")
    log:
        opj(OUTDIR, "logs/sratools/prefecth_{srr}.log")
    params:
        extra = "",
        accessions = "{srr}",
    threads: 8
    wrapper:
        get_wrapper('sra-tools', 'prefecth')

rule fasterq_dump:
    input:
        sra = opj(OUTDIR, "SRA/{srr}/{srr}.sra")
    output:
        expand(opj(OUTDIR, "raw/{{srr}}_{index}.fastq.gz"), index=[r[1] for r in RUN])
    log:
        opj(OUTDIR, "logs/sratools/fasterq-dump_{srr}.log")
    params:
        extra = "--split-3",
        # accession = "{srr}"
    threads: 8
    wrapper:
        get_wrapper('sra-tools', 'fasterq-dump')

rule fastq_merge:
    input:
        files = [opj(OUTDIR, "raw/{srr}_{index}.fastq.gz")]
    output:
        opj(OUTDIR, "raw/{srr}.fastq.gz")
    log:
        opj(OUTDIR, "logs/sratools/fastq-merge_{srr}.log")
    params:
        extra = ""
    threads: 8
    wrapper:
        get_wrapper('sra-tools', 'fastq-merge')