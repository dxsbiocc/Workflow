import json

rule PCA:
    input:
        exp = rules.merge_expression.output.exp
    output:
        pca = opj(OUTDIR, "expression/PCA.pdf"),
        qc = opj(OUTDIR, "expression/QC.pdf")
    log:
        opj(OUTDIR, "logs/diff/PCA.log")
    params:
        condition = json.dumps(GROUPS)
        height = 4,
        width = 6,
    threads: 5
    wrapper:
        get_wrapper('scripts', 'R', 'plot_PCA')

rule DEA:
    input:
        exp = rules.merge_expression.output.exp
    output:
        degs = opj(OUTDIR, "expression/DEGs.csv"),
        all_genes = opj(OUTDIR, "expression/DEGs.all.csv")
    log:
        opj(OUTDIR, "logs/diff/DEA.log")
    params:
        method = DIFF_TOOL,
        condition = json.dumps(GROUPS),
        pvalue = 0.05,
        log2FC = 1,
    threads: 5
    wrapper:
        get_wrapper('scripts', 'R', 'differential-expression')

rule volcano:
    input:
        all_genes = rules.DEA.output.all_genes
    output:
        volcano = opj(OUTDIR, "expression/volcano.pdf")
    log:
        opj(OUTDIR, "logs/diff/volcano.log")
    params:
        line_type = "line",
        height = 4,
        width = 6,
        color = "#2DB2EB #d8d8d8 #EB4232",
    wrapper:
        get_wrapper('scripts', 'R', 'plot_volcano')