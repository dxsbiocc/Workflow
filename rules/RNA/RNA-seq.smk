# configfile
configfile: os.path.join(PATH, "config/RNA-seq.yaml")

include: os.path.join(PATH, "rules/common/utils.smk")

READ_LENGTH = config['parameters']['star']['read_length']
# ------------------------- common rules ---------------------- #
# trimming
get_trimmed(TRIMMING)
# index
if MAPPING == 'star' and READ_LENGTH != 101:
    rule star_index:
        input:
            fasta = REF,
            gtf = GTF
        output:
            directory(INDEX),
        threads: 1
        params:
            extra = "",
            sjdbOverhang = READ_LENGTH - 1
        log:
            "logs/mapped/star_index.log",
        wrapper:
            get_wrapper("star", "index")
# mapping
get_mapping(MAPPING)
# dedup
get_dedup(DEDUP)
# ------------------------- special rules --------------------- #


# ---------------------------- outputs ------------------------ #
rule all:
    input:
        ""