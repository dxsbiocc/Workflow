import glob


if ASSEMBLE_MODE == 'solo':
    rule hifiasm_solo:
        input:
            seq = glob.glob('trimmed/*.filt.fastq.gz'),
        output:
            expand("assembly/genome.{suffix}.gfa", suffix=["p_ctg", "hap1.p_ctg", "hap2.p_ctg"]),
        log:
            "logs/hifiasm.log",
        params:
            extra = config['parameters']['hifiasm']['extra'],
            ultralong = "",
            primary = False,
        threads: config['parameters']['hifiasm']['threads']
        resources:
            mem_mb = config['parameters']['hifiasm']['resourses'],
        wrapper:
            get_wrapper("hifiasm")

elif ASSEMBLE_MODE == 'hic':
    rule hifiasm_hic:
        input:
            seq = glob.glob('trimmed/*.filt.fastq.gz'),
            hic1 = expand("trimmed/{hic}/{hic}.clean.R1.fq.gz", hic=DATA_DICT['HIC']),
            hic2 = expand("trimmed/{hic}/{hic}.clean.R2.fq.gz", hic=DATA_DICT['HIC']),
        output:
            expand("assembly/genome.{suffix}.gfa", suffix=["p_ctg", "hap1.p_ctg", "hap2.p_ctg"]),
        log:
            "logs/hifiasm_hic.log",
        params:
            extra = config['parameters']['hifiasm']['extra'],
            ultralong = "",
            primary = False,
        threads: config['parameters']['hifiasm']['threads']
        resources:
            mem_mb = config['parameters']['hifiasm']['resourses'],
        wrapper:
            get_wrapper("hifiasm")
elif ASSEMBLE_MODE == 'trio':
    rule yak_paternal:
        input:
            expand("trimmed/{trio}/{trio}.clean.R1.fq.gz")
        output:
            "yak/{trio}.yak"
        log:
            "logs/yak_{trio}.log"
        params:
            kmer = 8,
            extra = "-b37",
        threads: 8
        wrapper:
            get_wrapper("yak/count")

    rule hifiasm_trio:
        input:
            seq = glob.glob('trimmed/*.filt.fastq.gz'),
            yak = expand("yak/{trio}.yak", trio=TRIO),
        output:
            expand("assembly/genome.{suffix}.gfa", suffix=["p_ctg", "hap1.p_ctg", "hap2.p_ctg"]),
        log:
            "logs/hifiasm_trio.log",
        params:
            extra = config['parameters']['hifiasm']['extra'],
            ultralong = "",
            primary = False,
        threads: config['parameters']['hifiasm']['threads']
        resources:
            mem_mb = config['parameters']['hifiasm']['resourses'],
        wrapper:
            get_wrapper("hifiasm")
else:
    raise ValueError("`assembly_mode` must be one of 'solo', 'hic', 'trio'")