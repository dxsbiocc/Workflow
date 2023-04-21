rule genome_faidx:
    input:
        REFERENCE,
    output:
        REFERENCE + ".fai",
    log:
        "logs/genome/faidx.log",
    cache: True
    wrapper:
        get_wrapper('samtools', 'faidx')


rule genome_dict:
    input:
        REFERENCE,
    output:
        "resources/genome.dict",
    log:
        "logs/genome/create_dict.log",
    cache: True
    shell:
        get_wrapper('samtools', 'dict')


rule get_known_variation:
    input:
        # use fai to annotate contig lengths for GATK BQSR
        fai="resources/genome.fasta.fai",
    output:
        vcf="resources/variation.vcf.gz",
    log:
        "logs/check/get-known-variants.log",
    params:
        species = config["ref"]["species"],
        build = config["ref"]["build"],
        release = config["ref"]["release"],
        type = "all",
    cache: True
    wrapper:
        get_wrapper('reference', 'ensembl-variation')


rule remove_iupac_codes:
    input:
        "resources/variation.vcf.gz",
    output:
        "resources/variation.noiupac.vcf.gz",
    log:
        "logs/fix-iupac-alleles.log",
    conda:
        "../envs/rbt.yaml"
    cache: True
    shell:
        "rbt vcf-fix-iupac-alleles < {input} | bcftools view -Oz > {output}"


rule tabix_known_variants:
    input:
        "resources/variation.noiupac.vcf.gz",
    output:
        "resources/variation.noiupac.vcf.gz.tbi",
    log:
        "logs/check/tabix_variation.log",
    params:
        "-p vcf",
    cache: True
    wrapper:
        get_wrapper('tabix', 'index')

rule get_vep_cache:
    output:
        directory("resources/vep/cache"),
    params:
        species = config["ref"]["species"],
        build = config["ref"]["build"],
        release = config["ref"]["release"],
    log:
        "logs/vep/cache.log",
    wrapper:
        get_wrapper("vep", "cache")


rule get_vep_plugins:
    output:
        directory("resources/vep/plugins"),
    log:
        "logs/vep/plugins.log",
    params:
        release = config["ref"]["release"],
    wrapper:
        get_wrapper("vep", "plugins")
