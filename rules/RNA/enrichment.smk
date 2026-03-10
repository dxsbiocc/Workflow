rule enrichment:
    input:
        degs = rules.DEA.output.degs,
        db = PATHWAYDB,
    output:
        rda = opj(OUTDIR, "enrichment/enrichment.rda")
    log:
        opj(OUTDIR, "logs/enrichment/enrichment.log")
    params:
        species = config['parameters']['enrichment'].get('orgdb', 'org.Hs.eg.db'),
        go = config['parameters']['enrichment'].get('go'),
        pvalueCutoff = config['parameters']['enrichment'].get('pvalueCutoff', 0.05),
        qvalueCutoff = config['parameters']['enrichment'].get('qvalueCutoff', 0.2),
    wrapper:
        get_wrapper('scripts', 'R', 'enrichment')