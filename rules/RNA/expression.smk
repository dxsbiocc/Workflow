if DIFF_TOOL == "edger":
    rule edgeR:
        input:
            ""
        output:
            opj(OUTDIR, "")
        log:
            opj(OUTDIR, "logs/")
        threads: 5
        wrapper:
            get_wrapper('scripts', 'edger')
elif DIFF_TOOL == "deseq2":
    rule DESeq2:
        input:
            ""
        output:
            opj(OUTDIR, "")
        log:
            opj(OUTDIR, "")
        threads: 5
        wrapper:
            get_wrapper('scripts', 'deseq2')
elif DIFF_TOOL == "limma":
    rule limma:
        input:
            ""
        output:
            opj(OUTDIR, "")
        log:
            opj(OUTDIR, "")
        threads: 5
        wrapper:
            get_wrapper('scripts', 'limma')
else:
    raise ValueError(f'This analysis method: {DIFF_TOOL} not support!')

rule plot_DEGs:
    pass