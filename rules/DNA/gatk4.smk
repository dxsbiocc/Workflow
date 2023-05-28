rule baserecalibrator:
    input:
        bam = "dedup/{sample}/{sample}.rmdup.bam",
        ref = REFERENCE,
        dict = DICT,
        known = KNOWEN_SITE,  # optional known sites - single or a list
    output:
        recal_table = "variant/bqsr/recal_{sample}.grp",
    log:
        "logs/variant/gatk_haplotypecaller/{sample}.log",
    params:
        extra = "",  # optional
        java_opts = "",  # optional
    resources:
        mem_mb = 1024,
    wrapper:
        get_wrapper('gatk', 'baserecalibrator')

rule applybqsr:
    input:
        bam = "dedup/{sample}/{sample}.rmdup.bam",
        ref = REFERENCE,
        dict = DICT,
        recal_table = rules.baserecalibrator.output.recal_table,
    output:
        bam = "variant/bqsr/{sample}.bam",
    log:
        "logs/variant/gatk_applybqsr/{sample}.log",
    params:
        extra = "",  # optional
        java_opts = "",  # optional
        embed_ref = True,  # embed the reference in cram output
    resources:
        mem_mb = 1024,
    wrapper:
        get_wrapper('gatk', 'applybqsr')
        
rule call_variants:
    input:
        bam = rules.applybqsr.output.bam,
        ref = INDEX,
        idx = DICT,
        known = KNOWEN_SITE_DBSNP,
    output:
        gvcf = "variant/called/{sample}.{contig}.g.vcf.gz",
    log:
        "logs/variant/gatk_haplotypecaller/{sample}.{contig}.log",
    params:
        intervals = lambda wildcards: wildcards.contig,
        extra = "",
        java_opts = "",  # optional
    threads: 4
    resources:
        mem_mb = 1024,
    wrapper:
        get_wrapper('gatk', 'haplotypecaller')

rule combine_calls:
    input:
        ref = REFERENCE,
        gvcfs = expand(
            "variant/called/{sample}.{{contig}}.g.vcf.gz", sample=SAMPLES
        ),
    output:
        gvcf = "variant/called/all.{contig}.g.vcf.gz",
    log:
        "logs/gatk/combinegvcfs.{contig}.log",
    params:
        # intervals = lambda wildcards: wildcards.contig,
        extra = ""
    wrapper:
        get_wrapper('gatk', 'combinegvcfs')

rule genotypegvcfs:
    input:
        ref = REFERENCE,
        gvcf = rules.combine_calls.output.gvcf,
        # genomicsdb = rules.genomicsdbimport.output.db
    output:
        gvcf = "variant/genotyped/{contig}.g.vcf.gz",
    log:
        "logs/variant/gatk/genotypegvcfs/{contig}.log",
    wrapper:
        get_wrapper('gatk', 'genotypegvcfs')

rule merge_vcfs:
    input:
        vcfs = expand("variant/genotyped/{contig}.g.vcf.gz", contig=CONTIG),
    output:
        vcf = "variant/genotyped/all.vcf.gz",
    log:
        "logs/variant/picard/mergevcfs.log",
    params:
        extra = "",
    resources:
        mem_mb = 1024,
    wrapper:
        get_wrapper('picard', 'mergevcfs')

# select variant
rule selectvariants: 
    input: 
        vcf = rules.merge_vcfs.output,
        ref = REFERENCE
    output: 
        vcf = "variant/filtered/all.{mode}.raw.vcf.gz"
    log:
        "logs/gatk/selectvariants/{mode}.log",
    params:
        extra = lambda wildcards: "--select-type-to-include {}".format(wildcards.mode.upper()),
        java_opts = "",  # optional
    wrapper:
        get_wrapper('gatk', 'selectvariants')