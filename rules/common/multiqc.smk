rule multiqc_report:
    input:
        directory(OUTDIR)
    output:
        opj(OUTDIR, "report/report.html")
    params:
        extra = "",  # Optional: extra parameters for multiqc.
        use_input_files_only = True, # Optional, use only a.txt and don't search folder samtools_stats for files
    log:
        opj(OUTDIR, "logs/report/multiqc.log")
    wrapper:
        get_wrapper('multiqc')