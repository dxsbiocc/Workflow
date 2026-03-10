if REMOVE_HOST:
    rule kneaddata:
        input:
            fq_list = get_fastq,
            db = config['parameters']['host']['db'],
        output:
            fastq = opj(OUTDIR, "kneaddata/{sample}/{sample}_kneaddata.fastq"),
            contam = opj(OUTDIR, "kneaddata/{sample}/{sample}_contam.fastq"),
        log:
            opj(OUTDIR, "logs/kneaddata/{sample}.kneaddata.log"),
        params:
            extra = config['parameters']['host']['extra'],
        threads: config['parameters']['host']['threads'],
        wrapper:
            get_wrapper('kneaddata')
else:
    rule concat_reads:
        input:
            fastq = expand(opj(OUTDIR, "trimmed/{{sample}}/{{sample}}.clean.{run}.fq.gz"), run=RUN),
        output:
            opj(OUTDIR, "concat_fastq/{sample}.merged.fq.gz"),
        log:
            opj(OUTDIR, "logs/concat_fastq/{sample}.merged.log"),
        threads: 8
        shell:
            "cat {input} > {output}"

if READS and isinstance(READS, int):
    rule downsample_reads:
        input:
            fastq = rules.kneaddata.output.fastq if REMOVE_HOST else rules.concat_reads.output,
        output:
            opj(OUTDIR, "downsample/{sample}.sampled.fq.gz"),
        log:
            opj(OUTDIR, "logs/downsample/{sample}.sampled.log"),
        params:
            n = READS,
            seed = 12345,
        wrapper:
            get_wrapper("seqtk", "subsample")
else:
    rule all_reads:
        input:
            fastq = rules.kneaddata.output.fastq if REMOVE_HOST else rules.concat_reads.output,
        output:
            opj(OUTDIR, "downsample/{sample}.sampled.fq.gz")
        log:
            opj(OUTDIR, "logs/downsample/{sample}.sampled.log"),
        threads: 8
        shell:
            "ln -s {input} {output}"