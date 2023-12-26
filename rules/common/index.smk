if MAPPING == 'bwa':
    rule bwa_index:
        input:
            fasta = REFERENCE
        output:
            idx = multiext(
                INDEX, 
                ".amb", 
                ".ann", 
                ".bwt", 
                ".pac", 
                ".sa"
            ),
        threads: 20
        params:
            extra = ""
        log:
            opj(OUTDIR, "logs/mapped/bwa_index.log"),
        wrapper:
            get_wrapper("bwa", "index")
elif MAPPING == 'bwa-mem2':
    rule bwa_mem2_index:
        input:
            fasta = REFERENCE
        output:
            multiext(
                INDEX, 
                ".0123", 
                ".amb", 
                ".ann", 
                ".pac", 
                ".bwt.2bit.64"
            ),
        threads: 20
        params:
            extra = ""
        log:
            opj(OUTDIR, "logs/mapped/bwa_index.log"),
        wrapper:
            get_wrapper("bwa-mem2", "index")
elif MAPPING == 'bowtie2':
    rule bowtie2_index:
        input:
            fasta = REFERENCE
        output:
            multiext(
                INDEX,
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            ),
        threads: 20
        params:
            extra = ""
        log:
            opj(OUTDIR, "logs/mapped/bowtie2_index.log"),
        wrapper:
            get_wrapper("bowtie2", "index")
elif MAPPING == 'hisat2':
    rule hisat2_index:
        input:
            fasta = REFERENCE
        output:
            directory(INDEX),
        threads: 20
        params:
            extra = "",
            prefix = GENOME
        log:
            opj(OUTDIR, "logs/mapped/hisat2_index.log"),
        wrapper:
            get_wrapper("hisat2", "index")
elif MAPPING == 'star':
    rule star_index:
        input:
            fasta = REFERENCE,
            gtf = GTF
        output:
            directory(INDEX),
        threads: 20
        params:
            extra = "",
            sjdbOverhang = config['parameters']['star']['read_length'] - 1
        log:
            opj(OUTDIR, "logs/mapped/star_index.log"),
        wrapper:
            get_wrapper("star", "index")
elif MAPPING == 'minimap2':
    rule minimap2_index:
        input:
            fasta = REFERENCE
        output:
            INDEX,
        threads: 20
        params:
            extra = ""
        log:
            opj(OUTDIR, "logs/mapped/minimap2_index.log"),
        wrapper:
            get_wrapper("minimap2", "index")
else:
    raise ValueError("Unknown mapping tool: {}".format(MAPPING))