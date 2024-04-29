sink(snakemake@log[[1]])
options(warn=-1)
# library packages
suppressMessages(library(tidyverse, quietly = TRUE))
suppressMessages(library(clusterProfiler, quietly = TRUE))

enrichment <- function(degs, db = NULL, go = NULL, species = NULL, 
                       pvalueCutoff = 0.05, qvalueCutoff = 0.2) {
    DEGs <- read_csv(degs, show_col_types = FALSE)
    DEGs.Up <- dplyr::filter(DEGs, direction == "Up")
    DEGs.Down <- dplyr::filter(DEGs, direction == "Down")
    result <- list()
    if (!is.null(go)) {
        if (is.null(species)) stop("species must be set.")
        if (nrow(DEGs.Up) > 0) {
            cat("GO Enrichment for upregalated gene!", sep = "\n")
            result$up.ora.go <- enrichmentGO(DEGs.Up, go, species, pvalueCutoff, qvalueCutoff)
            result$up.gsea.go <- enrichmentGO(DEGs.Up, go, species, pvalueCutoff, qvalueCutoff, TRUE)
        }
        if (nrow(DEGs.Down) > 0) {
            cat("GO Enrichment for downregalated gene!", sep = "\n")
            result$down.ora.go <- enrichmentGO(DEGs.Down, go, species, pvalueCutoff, qvalueCutoff)
            result$down.gsea.go <- enrichmentGO(DEGs.Down, go, species, pvalueCutoff, qvalueCutoff, TRUE)
        }
    }
    if (!is.null(db)) {
        if (nrow(DEGs.Up) > 0) {
            cat("DB Enrichment for upregalated gene!", sep = "\n")
            result$up.ora.db <- enrichmentDB(DEGs.Up, db, pvalueCutoff, qvalueCutoff)
            result$up.gsea.db <- enrichmentDB(DEGs.Up, db, pvalueCutoff, qvalueCutoff, TRUE)
        }
        if (nrow(DEGs.Down) > 0) {
            cat("DB Enrichment for downregalated gene!", sep = "\n")
            result$down.ora.db <- enrichmentDB(DEGs.Down, db, pvalueCutoff, qvalueCutoff)
            result$down.gsea.db <- enrichmentDB(DEGs.Down, db, pvalueCutoff, qvalueCutoff, TRUE)
        }
    }
    result
}

enrichmentGO <- function(DEGs, go, species, pvalueCutoff, qvalueCutoff, gsea = FALSE) {
    if (dim(DEGs)[1] < 10) {
        return(tibble())
    }
    gene.id <- suppressMessages(bitr(
        DEGs$symbol, fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = species
    )) %>%
        mutate(ENTREZID = as.integer(ENTREZID)) %>%
        inner_join(DEGs, by = c("SYMBOL" = "symbol"))
    if (gsea) {
        geneList <- sort(structure(gene.id$log2FC, names = gene.id$ENTREZID), decreasing = TRUE)
        res <- suppressMessages(gseGO(geneList, OrgDb = species, ont = go, keyType = "ENTREZID", 
                                      pvalueCutoff = pvalueCutoff))
    } else {
        res <- suppressMessages(enrichGO(gene = gene.id$ENTREZID, OrgDb = species, ont = go, 
                                         pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff))
    }
    return(res)
}

enrichmentDB <- function(DEGs, db, pvalueCutoff, qvalueCutoff, gsea = FALSE) {
    if (endsWith(db, "csv")) {
        TERM2GENE <- read_csv(db)
    } else if (endsWith(db, "gmt")) {
        TERM2GENE <- read.gmt(db)
    }
    else {
        stop("unsupported pathway file")
    }
    if (gsea) {
        geneList <- sort(structure(DEGs$log2FC, names = DEGs$symbol), decreasing = TRUE)
        res <- suppressMessages(GSEA(geneList, TERM2GENE = TERM2GENE, pvalueCutoff = pvalueCutoff))
    } else {
        res <- suppressMessages(enricher(gene = DEGs$symbol, TERM2GENE = TERM2GENE,
                                         pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff))
    }
    return(res)
}

# parse arguments
degs <- snakemake@input[["degs"]]
db <- snakemake@input[["db"]]
species <- as.character(snakemake@params[["species"]])
go <- snakemake@params[["go"]]
outfile <- snakemake@output[['rda']]
pvalueCutoff <- as.numeric(snakemake@params[["pvalue"]])
qvalueCutoff <- as.numeric(snakemake@params[["qvalue"]])

# library
cat(paste0("Library package: ", species))
suppressMessages(library(species, quietly = TRUE, character.only = TRUE))
# runing
res <- enrichment(degs, db, go, species, pvalueCutoff, qvalueCutoff)
cat(paste0("Save result to ", outfile, "."), sep = "\n")
save(res, file = outfile)
# end
sink()
