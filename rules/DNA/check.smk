import os

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

if not VEP_CACHE:
    VEP_CACHE = "resources/vep/cache"
    rule get_vep_cache:
        output:
            directory("resources/vep/cache"),
        params:
            species = VEP_SPECIES,
            build = VEP_BUILD,
            release = VEP_RELEASE,
        log:
            "logs/vep/cache.log",
        wrapper:
            get_wrapper("vep", "cache")

if not VEP_PLUGINS:
    VEP_PLUGINS = "resources/vep/plugins"
    rule get_vep_plugins:
        output:
            directory("resources/vep/plugins"),
        log:
            "logs/vep/plugins.log",
        params:
            release = VEP_RELEASE,
        wrapper:
            get_wrapper("vep", "plugins")
