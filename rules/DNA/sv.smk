rule delly_call:
    input:
        alns = "dedup/{sample}/{sample}.rmdup.bam",
        ref = REFERENCE,
        # optional
        exclude = DELLY_EXCLUDE,
    output:
        "sv/{sample}.calls.vcf.gz",
    params:
        extra = "",  # optional parameters for delly (except -g, -x)
    log:
        "logs/sv/delly/call_{sample}.log",
    threads: 2  # It is best to use as many threads as samples
    wrapper:
        get_wrapper('delly')
        
# merge variant    
rule merge_sv:
    input:
        vcfs = expand("sv/{sample}.calls.vcf.gz", sample=SAMPLES),
    output:
        vcf = "sv/all.vcf.gz",
    log:
        "logs/sv/picard/merge-sv.log",
    wrapper:
        get_wrapper('picard', 'mergevcfs')