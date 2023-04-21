def get_file_list():
    return []


rule multiqc_report:
    input:
        get_file_list
    output:
        "report/report.html"
    params:
        extra = "",  # Optional: extra parameters for multiqc.
        use_input_files_only = True, # Optional, use only a.txt and don't search folder samtools_stats for files
    log:
        "logs/report/multiqc.log"
    wrapper:
        get_wrapper('multiqc')