if MAPPING == 'bwa':
    rule bwa_index:
        input:
            fasta = REF
        output:
            idx = multiext(
                os.path.splitext(REF)[0], 
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
            "logs/mapped/bwa_index.log",
        wrapper:
            get_wrapper("bwa", "index")
elif MAPPING == 'bwa-mem2':
    rule bwa_mem2_index:
        input:
            fasta = REF
        output:
            multiext(
                os.path.splitext(REF)[0], 
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
            "logs/mapped/bwa_index.log",
        wrapper:
            get_wrapper("bwa-mem2", "index")
elif MAPPING == 'bowtie2':
    rule bowtie2_index:
        input:
            fasta = REF
        output:
            multiext(
                os.path.splitext(REF)[0],
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
            "logs/mapped/bowtie2_index.log",
        wrapper:
            get_wrapper("bowtie2", "index")
elif MAPPING == 'hisat2':
    rule hisat2_index:
        input:
            fasta = REF
        output:
            directory(os.path.dirname(REF)),
        threads: 20
        params:
            extra = ""
        log:
            "logs/mapped/hisat2_index.log",
        wrapper:
            get_wrapper("hisat2", "index")
elif MAPPING == 'star':
    rule star_index:
        input:
            fasta = REF,
            gtf = GTF
        output:
            directory(os.path.dirname(REF)),
        threads: 20
        params:
            extra = "",
            sjdbOverhang = 100
        log:
            "logs/mapped/star_index.log",
        wrapper:
            get_wrapper("star", "index")
elif MAPPING == 'minimap2':
    rule minimap2_index:
        input:
            fasta = REF
        output:
            f"{os.path.splitext(REF)[0]}.mmi",
        threads: 20
        params:
            extra = ""
        log:
            "logs/mapped/minimap2_index.log",
        wrapper:
            get_wrapper("minimap2", "index")
else:
    raise ValueError("Unknown mapping tool: {}".format(MAPPING))