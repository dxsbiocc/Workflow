rule vep_annotate:
    input:
        calls = rules.merge_calls.output.vcf,
        cache = VEP_CACHE,
        plugins = VEP_PLUGINS_LOCAL,
    output:
        calls = report(
            opj(OUTDIR, "annotated/all.vcf.gz"),
            caption = "report/vcf.rst",
            category = "Calls",
        ),
        stats = report(
            opj(OUTDIR, "annotated/stats/all.stats.html"),
            caption = "report/stats.rst",
            category = "Calls",
        ),
    params:
        plugins = VEP_PLUGINS,
        extra = "",
    log:
        opj(OUTDIR, "logs/vep/annotate.log"),
    threads: 4
    wrapper:
        get_wrapper('vep', 'annotate')