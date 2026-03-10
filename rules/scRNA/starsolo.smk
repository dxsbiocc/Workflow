if PLATFORM == '10x':
    rule starsolo_10x:
        input:
            fq1 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R1.fq.gz"),
            fq2 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R2.fq.gz"),
            index = INDEX,
        output:
            directory(opj(OUTDIR, "starsolo/{sample}"))
        log:
            opj(OUTDIR, "logs/mapped/STARsolo_{sample}.log"),
        params:
            # optional parameters
            scrna = PLATFORM,
            barcode = os.path.join(DATABASE, 'barcode'),
            extra = config['parameters']['star']['extra'] + \
                    " --outSAMattrRGline ID:{sample} SM:{sample}",
        threads: config['parameters']['star']['threads']
        wrapper:
            get_wrapper("star", "solo")
elif PLATFORM == 'drop':
    rule starsolo_drop:
        input:
            fq1 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R1.fq.gz"),
            fq2 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R2.fq.gz"),
            index = INDEX,
        output:
            directory(opj(OUTDIR, "starsolo/{sample}"))
        log:
            opj(OUTDIR, "logs/mapped/STARsolo_{sample}.log"),
        params:
            # optional parameters
            scrna = PLATFORM,
            barcode = os.path.join(DATABASE, 'barcode'),
            cblen = 12,
            umilen = 8,
            extra = config['parameters']['star']['extra'] + \
                    " --outSAMattrRGline ID:{sample} SM:{sample}",
        threads: config['parameters']['star']['threads']
        wrapper:
            get_wrapper("star", "solo")
elif PLATFORM == 'indrop':
    rule starsolo_indrop:
        input:
            fq1 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R1.fq.gz"),
            fq2 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R2.fq.gz"),
            index = INDEX,
        output:
            directory(opj(OUTDIR, "starsolo/{sample}"))
        log:
            opj(OUTDIR, "logs/mapped/STARsolo_{sample}.log"),
        params:
            # optional parameters
            scrna = PLATFORM,
            barcode = os.path.join(DATABASE, 'barcode'),
            adapter = 'GAGTGATTGCTTGTGACGCCTT',                   # these could be GAGTGATTGCTTGTGACGCCTT or GAGTGATTGCTTGTGACGCCAA 
            extra = config['parameters']['star']['extra'] + \
                    " --outSAMattrRGline ID:{sample} SM:{sample}"
        threads: config['parameters']['star']['threads']
        wrapper:
            get_wrapper("star", "solo")
elif PLATFORM == 'smart2':
    rule starsolo_smart2:
        input:
            fq1 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R1.fq.gz"),
            fq2 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R2.fq.gz"),
            index = INDEX,
        output:
            directory(opj(OUTDIR, "starsolo/{sample}"))
        log:
            opj(OUTDIR, "logs/mapped/STARsolo_{sample}.log"),
        params:
            # optional parameters
            scrna = PLATFORM,
            extra = config['parameters']['star']['extra'] + \
                    " --outSAMattrRGline ID:{sample} SM:{sample}"
        threads: config['parameters']['star']['threads']
        wrapper:
            get_wrapper("star", "solo")
elif PLATFORM == 'cel':
    rule starsolo_cel:
        input:
            fq1 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R1.fq.gz"),
            fq2 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R2.fq.gz"),
            index = INDEX,
        output:
            directory(opj(OUTDIR, "solo/{sample}")),
        log:
            opj(OUTDIR, "logs/mapped/STARsolo_{sample}.log"),
        params:
            # optional parameters
            scrna = PLATFORM,
            cblen = 6,
            umilen = 6,
            extra = config['parameters']['star']['extra'] + \
                    " --outSAMattrRGline ID:{sample} SM:{sample}"
        threads: config['parameters']['star']['threads']
        wrapper:
            get_wrapper("star", "solo")
elif PLATFORM == 'custom':
    rule starsolo_custom:
        input:
            fq1 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R1.fq.gz"),
            fq2 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R2.fq.gz"),
            index = INDEX,
        output:
            directory(opj(OUTDIR, "solo/{sample}")),
        log:
            opj(OUTDIR, "logs/mapped/STARsolo_{sample}.log"),
        params:
            # optional parameters
            scrna = PLATFORM,
            extra = config['parameters']['star']['extra'] + \
                    " --outSAMattrRGline ID:{sample} SM:{sample}"
        threads: config['parameters']['star']['threads']
        wrapper:
            get_wrapper("star", "solo")
else:
    raise ValueError(f"Unsupported PLATFORM: {PLATFORM}")