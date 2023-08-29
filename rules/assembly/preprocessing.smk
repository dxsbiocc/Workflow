# HIFI adapter filter
rule hifiadapterfit:
    input:
        seq = lambda wildcards: DATA.loc[wildcards.hifi, "fastq1"],
        db = get_adapter('HIFI'),
        script = get_script('pbadapterfilt.sh')
    output:
        directory("trimmed/{hifi}")
    params:
        min_length = 44,
        min_percentage = 97,
    threads: 10
    log:
        "logs/hifiadapterfit/{hifi}.log"
    wrapper:
        get_wrapper("hifiadapterfit")


# NGS reads quality control
rule fastp:
    input:
        reads = get_fastq,
    output:
        trimmed = expand("trimmed/{{sample}}/{{sample}}.clean.{run}.fq.gz", run=RUN),
        # or in a single file
        # unpaired = "trimmed/{sample}/{sample}.singletons.fastq",
        # merged = "fastp/{sample}.merged.fastq",
        failed = "trimmed/{sample}/{sample}.failed.fastq",
        html = "trimmed/report/{sample}.fastp.html",
        json = "trimmed/report/{sample}.fastp.json"
    log:
        "logs/trimmed/fastp_{sample}.log"
    params:
        extra = "-g -q 30 -u 50 -n 15 -l 36 -w 4",  # optional parameters
        adapters = "--adapter_fasta " + os.path.join(PATH, 'data/adapters/NexteraPE-PE.fa')
    threads: 2
    wrapper:
        get_wrapper('fastp')