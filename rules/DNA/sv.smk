rule delly_call:
    input:
        alns = opj(OUTDIR, "dedup/{sample}/{sample}.rmdup.bam"),
        ref = REFERENCE,
        # optional
        exclude = DELLY_EXCLUDE,
    output:
        opj(OUTDIR, "sv/{sample}.calls.vcf.gz"),
    params:
        extra = "",  # optional parameters for delly (except -g, -x)
    log:
        opj(OUTDIR, "logs/sv/delly/call_{sample}.log"),
    threads: 2  # It is best to use as many threads as samples
    wrapper:
        get_wrapper('delly')
        
# merge variant    
rule merge_sv:
    input:
        vcfs = expand(opj(OUTDIR, "sv/{sample}.calls.vcf.gz"), sample=SAMPLES),
    output:
        vcf = opj(OUTDIR, "sv/all.vcf.gz"),
    log:
        opj(OUTDIR, "logs/sv/picard/merge-sv.log"),
    wrapper:
        get_wrapper('picard', 'mergevcfs')