options(warn=-1)
# output to log file
sink(snakemake@log[[1]])

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(jsonlite)))
suppressWarnings(suppressMessages(library(ggrepel)))

# QC
plot_QC <- function(exprs, outfile, height, width) {
    par(cex = 0.7)
    n.sample <- ncol(exprs)
    if(n.sample > 40) par(cex = 0.5)
    cols <- rainbow(n.sample * 1.2)
    pdf(outfile, height = height, width = width)
    boxplot(log2(exprs + 1), col = cols, main = "log2(FPKM + 1)", las = 2)
    dev.off()
}

# PCA
plot_PCA <- function(exprs, coldata, outfile, height, width) {
    pca_result <- prcomp(exprs, center = TRUE, scale. = TRUE)
    pca.data <- data.frame(sample = rownames(pca_result$rotation),
                        group = coldata[colnames(exprs), "group"],
                        pca_result$rotation)

    p <- ggplot(pca.data, aes(x = PC1, y = PC2)) +
        geom_point(aes(fill = group), shape = 21, size = 5) +
        geom_text_repel(aes(label = sample)) +
        labs(title = "PCA Plot", x = "PC1", y = "PC2") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    pdf(outfile, height = height, width = width)
    print(p)
    dev.off()
}

main <- function(infile, condition, outpca, logfile, outqc = NULL, height, width) {
    writeLines("parsing condition ...", logfile)
    condition <- fromJSON(condition)
    coldata <- do.call(rbind, lapply(names(condition), function(x) {
        data.frame(
            sample = condition[[x]],
            group = x
        )
    })) %>%
        mutate(group = factor(group)) %>%
        tibble::column_to_rownames("sample")
    cat(paste("read FPKM from", infile), sep = "\n")
    fpkm <- suppressMessages(read_tsv(infile)) %>%
        dplyr::select(1, ends_with("FPKM")) %>%
        tibble::column_to_rownames("gene_name") %>%
        rename_with(function(x) gsub("_FPKM", "", x))

    if (!is.null(outqc)) {
        cat(paste("plot QC boxplot", outqc), sep = "\n")
        plot_QC(fpkm, outqc, height, width)
    }
    cat(paste("plot PCA to", outpca), sep = "\n")
    plot_PCA(fpkm, coldata, outpca, height, width)
    cat("PCA done!", sep = "\n")
}

# parse arguments
infile <- snakemake@input[["exp"]]
outpca <- snakemake@output[["pca"]]
outqc <- snakemake@output[["qc"]]
logfile <- snakemake@log[[1]]
condition <- snakemake@params[["condition"]]

height <- as.numeric(snakemake@params[["height"]])
width <- as.numeric(snakemake@params[["width"]])

# main
main(infile, condition, outpca, logfile, outqc, height, width)

sink()