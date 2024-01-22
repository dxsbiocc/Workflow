if FUSION == "arriba":
    rule arriba:
        input:
            # STAR bam containing chimeric alignments
            bam = opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam"),
            # path to reference genome
            genome = REFERENCE,
            # path to annotation gtf
            annotation = GTF,
        output:
            # approved gene fusions
            fusions = opj(OUTDIR, "fusion/arriba/{sample}.tsv"),
            # discarded gene fusions
            discarded = opj(OUTDIR, "fusion/arriba/{sample}.discarded.tsv") # optional
        log:
            opj(OUTDIR, "logs/fusion/arriba_{sample}.log")
        threads: 10
        params:
            # arriba blacklist file
            blacklist = BLACKLIST, # strongly recommended
            # file containing known fusions
            known_fusions = KNOWN_FUSIONS, # optional
            # file containing information from structural variant analysis
            sv_file = SV_FILE, # optional
            # optional parameters
            extra = config['parameters']['arriba']['extra']
        wrapper:
            get_wrapper("arriba")
elif FUSION == "star-fusion":
    if MAPPING == "star":
        rule star_fusion_kickstart:
            input:
                junction = opj(OUTDIR, "mapped/{sample}/Chimeric.out.junction"),
                index = GENOME_LIB,
            output:
                directory(opj(OUTDIR, "fusion/star-fusion/{sample}"))
            log:
                opj(OUTDIR, "logs/fusion/star_fusion_kickstart_{sample}.log")
            threads: 10
            params:
                extra = config['parameters']['star-fusion']['extra']
            wrapper:
                get_wrapper("star-fusion")
    else:
        rule star_fusion:
            input:
                reads = expand(opj(OUTDIR, "trimmed/{{sample}}/{{sample}}.clean.{run}.fq.gz"), run=RUN),
                index = GENOME_LIB,
            output:
                directory(opj(OUTDIR, "fusion/star-fusion/{sample}"))
            log:
                opj(OUTDIR, "logs/fusion/star_fusion_{sample}.log")
            threads: 8
            params:
                extra = config['parameters']['star-fusion']['extra']
            wrapper:
                get_wrapper("star-fusion")
else:
    # raise ValueError(f"{FUSION} not support! Only used `arriba` or `star-fusion`!")
    pass