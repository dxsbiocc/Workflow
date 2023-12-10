if MAPPING == "star":
    rule star_fusion_kickstart:
        input:
            junction = opj(OUTDIR, "mapped/{sample}/Chimeric.out.junction"),
            index = GENOME_LIB,
        output:
            directory(opj(OUTDIR, "fusion/{sample}/"))
        log:
            opj(OUTDIR, "logs/fusion/star_fusion_kickstart_{sample}.log")
        threads: 10
        wrapper:
            get_wrapper("star-fusion")
else:
    rule star_fusion:
        input:
            reads = expand(opj(OUTDIR, "trimmed/{{sample}}/{{sample}}.clean.{run}.fq.gz"), run=RUN),
            index = GENOME_LIB,
        output:
            directory(opj(OUTDIR, "fusion/{sample}/"))
        log:
            opj(OUTDIR, "logs/fusion/star_fusion_{sample}.log")
        threads: 8
        wrapper:
            get_wrapper("star-fusion")