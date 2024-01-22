sink(snakemake@log[[1]])
options(warn=-1)
# library packages
suppressWarnings(suppressMessages(library("edgeR")))
suppressWarnings(suppressMessages(library("limma")))
suppressWarnings(suppressMessages(library("jsonlite")))
suppressWarnings(suppressMessages(library("tidyverse")))
suppressWarnings(suppressMessages(library("DESeq2")))

parse_condition <- function(condition) {
    condition <- '{"Cd": ["Cd-1", "Cd-2", "Cd-3"], "Ctrl": ["Ctrl-1", "Ctrl-2", "Ctrl-3"]}'
    condition <- fromJSON(condition)

    ctrl <- c("Ctrl", "ctrl", "Control", "control", "Normal", "normal")
    groups <- names(condition)
    index <- groups %in% ctrl
    level <- c(groups[index], groups[!index])

    coldata <- do.call(rbind, lapply(names(condition), function(x) {
        data.frame(
            sample = condition[[x]],
            group = x
        )
    })) %>%
        mutate(group = factor(group, levels = level)) %>%
        tibble::column_to_rownames("sample")

    return(coldata)
}

dea_deseq2 <- function(counts_matrix, coldata) {
    dds <- DESeqDataSetFromMatrix(
        countData = counts_matrix,
        colData = coldata,
        design = ~group
    )
    dds <- dds[rowSums(counts(dds)) > 0, ]
    dds <- DESeq(dds)
    res <- results(dds)

    plot_data <- as.data.frame(res) %>%
        tibble::rownames_to_column("symbol") %>%
        dplyr::select(c(1, 3, 7)) %>%
        dplyr::rename_with(~ c("symbol", "log2FC", "p.value"))
    return(plot_data)
}

dea_limma <- function(counts_matrix, coldata) {
    groups <- coldata[colnames(counts_matrix), ]
    nf <- calcNormFactors(counts_matrix)
    design <- model.matrix(~groups)
    y <- voom(counts_matrix, design, lib.size = colSums(counts_matrix) * nf)
    fit <- lmFit(y, design)
    fit <- eBayes(fit)
    res <- topTable(fit, coef = 2, n = length(counts_matrix[, 1]), sort = "p")

    plot_data <- as.data.frame(res) %>%
        tibble::rownames_to_column("symbol") %>%
        dplyr::select(c(1, 2, 6)) %>%
        dplyr::rename_with(~ c("symbol", "log2FC", "p.value"))
    return(plot_data)
}

dea_edger <- function(counts_matrix, coldata, method = "glm") {
    groups <- coldata[colnames(counts_matrix), ]
    dgelist <- DGEList(counts = counts_matrix, group = groups)
    # filter low count, using CPM
    keep <- rowSums(cpm(dgelist) > 1) >= 2
    dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
    # TMM normalization
    dgelist_norm <- calcNormFactors(dgelist, method = "TMM")
    # grouping
    design <- model.matrix(~groups)
    # Estimating the dispersion of gene expression values
    dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
    # differentail expression analysis
    if (method == "glm") {
        # 1. Negative binomial generalized log-linear models
        fit <- glmFit(dge, design, robust = TRUE)
        lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
    } else if (method == "glmQL") {
        # 2. Proposed Likelihood Negative Binomial Generalized Log Linear Models
        fit <- glmQLFit(dge, design, robust = TRUE)
        lrt <- topTags(glmQLFTest(fit), n = nrow(dgelist$counts))
    } else {
        return(NULL)
    }

    plot_data <- as.data.frame(lrt$table) %>%
        tibble::rownames_to_column("symbol") %>%
        dplyr::select(c(1, 2, 6)) %>%
        dplyr::rename_with(~ c("symbol", "log2FC", "p.value"))
    return(plot_data)
}

dea <- function(infile, condition, outfile, out_all, logfile, methods = "deseq2", pvalue = 0.05, log2FC = 1) {
    coldata <- parse_condition(condition)

    counts_matrix <- suppressMessages(read_tsv(infile)) %>%
        dplyr::select(1, ends_with("counts")) %>%
        tibble::column_to_rownames("gene_name") %>%
        mutate_all(as.integer) %>%
        rename_with(function(x) gsub("_counts", "", x))

    cat(paste("running:", methods), sep = "\n")
    if (methods == "deseq2") {
        result <- suppressMessages(dea_deseq2(counts_matrix, coldata))
    } else if (methods == "edger") {
        result <- suppressMessages(dea_edger(counts_matrix, coldata, method = "glm"))
    } else if (methods == "limma") {
        result <- suppressMessages(dea_limma(counts_matrix, coldata))
    } else {
        cat(paste("unknown method:", methods), sep = "\n")
        stop()
    }
    cat("Differential expressional analysis done!", sep = "\n")
    # write output
    cat(paste("All gene result writing to", out_all), sep = "\n")
    write.csv(result, file = out_all, row.names = FALSE)

    dges <- dplyr::filter(result, abs(log2FC) > 1 & p.value < pvalue) %>%
        mutate(direction = ifelse(log2FC > 1, "Up", "Down"))
    cat(paste("Differential expressional genes writing to", outfile), sep = "\n")
    write.csv(dges, file = outfile, row.names = FALSE)
}

# parse arguments
infile <- snakemake@input[["exp"]]
out_degs <- snakemake@output[["degs"]]
out_all <- snakemake@output[["all_genes"]]
logfile <- snakemake@log[[1]]
condition <- snakemake@params[["condition"]]
method <- snakemake@params[["method"]]
pvalue <- as.numeric(snakemake@params[["pvalue"]])
log2FC <- as.numeric(snakemake@params[["log2FC"]])

# runing
dea(infile, condition, out_degs, out_all, logfile, method, pvalue, log2FC)
# end
sink()
