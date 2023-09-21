# filter variant, hard and VQSR filter
if FILTER == 'hard':  # hard filter
    rule hard_filters:
        input:
            vcf = rules.selectvariants.output.vcf,
            ref = REFERENCE
        output:
            vcf = "variant/filtered/all.{mode}.filtered.vcf.gz",
        params:
            filters = lambda wildcard: get_filter(wildcard.mode), 
        log:
            "logs/gatk/variantfiltration/{mode}.log",
        resources:
            mem_mb=1024,
        wrapper:
            get_wrapper('gatk', 'variantfiltration')
            
else: # VQSR filter
    rule variantrecalibrator:
        input:
            vcf = rules.selectvariants.output.vcf,
            ref = REFERENCE,
            dict = DICT,
            mills = MILLS,
            mills_idx = MILLS_IDX,
            omni = OMNI,
            omni_idx = OMIN_IDX,
            g1k = G1K,
            g1k_idx = G1K_IDX,
            dbsnp = DBSNP,
            dbsnp_idx = DBSNP_IDX
        output:
            recal = "variant/filtered/all.{mode}.vqsr.recal",
            tranches = "variant/filtered/all.{mode}.tranches",
        params:
            mode = lambda wildcards: '{}'.format(wildcards.mode).upper(),
            resources={
                "mills": {"known": False, "training": True, "truth": True, "prior": 15.0},
                "omni": {"known": False, "training": True, "truth": False, "prior": 12.0},
                "g1k": {"known": False, "training": True, "truth": False, "prior": 10.0},
                "dbsnp": {"known": True, "training": False, "truth": False, "prior": 2.0},
            },
            annotation=["MQ", "QD", "SB"],
            extra="--max-gaussians 2",  # optional
        threads: 1
        resources:
            mem_mb = 1024,
        log:
            "logs/gatk/variantrecalibrator/{mode}.log",
        wrapper:
            get_wrapper('gatk', 'variantrecalibrator')

    rule apply_vqsr:
        input:
            vcf = rules.selectvariants.output,
            recal = rules.variantrecalibrator.output.vcf,
            tranches = rules.variantrecalibrator.output.tranches,
            ref = REFERENCE,
        output:
            vcf = "variant/filtered/all.{mode}.filtered.vcf.gz",
        params:
            mode = lambda wildcards: '{}'.format(wildcards.mode).upper(),
            extra = "",  # optional
        resources:
            mem_mb = 50,
        log:
            "logs/gatk/applyvqsr/{mode}.log",
        wrapper:
            get_wrapper('gatk', 'applyvqsr')
            
# merge variant    
rule merge_calls:
    input:
        vcfs = expand("variant/filtered/all.{mode}.filtered.vcf.gz", mode=MODE),
    output:
        vcf = "variant/filtered/all.vcf.gz",
    log:
        "logs/picard/merge-filtered.log",
    wrapper:
        get_wrapper('picard', 'mergevcfs')