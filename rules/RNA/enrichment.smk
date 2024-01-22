rule enrichment_go:
    input:
        "path/to/diff_genes.txt"
    output:
        "path/to/enrichment_results.txt"
    script:
        """
        library(clusterProfiler)
        diff_genes <- read.table({{input}}, header = TRUE)
        ego <- enrichGO(gene = diff_genes$gene_id,
                        universe = universe,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENSEMBL",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05,
                        readable = TRUE)
        write.table(ego, file = {{output}}, sep = "\t")
        """

rule enrichment_kegg:
    pass

rule enrichment_MSigDB:
    pass

rule gsea_go:
    pass

rule gsea_kegg:
    pass

rule gsea_MSigDB:
    pass