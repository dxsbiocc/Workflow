rule haplotypecaller:
    input:
        bam = get_sample_bams,
        ref = REFERENCE,
        idx = "resources/genome.dict",
        known = KNOWN_SNP[0],
    output:
        gvcf = protected("called/{sample}.g.vcf.gz"),
    log:
        "logs/called/haplotypecaller_{sample}.log",
    params:
        extra = get_call_variants_params,
        java_opts = "",  # optional
    wrapper:
        get_wrapper("gatk", "haplotypecaller")