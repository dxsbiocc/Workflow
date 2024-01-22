sink(snakemake@log[[1]])
options(warn = -1)

suppressMessages(library(tidyverse))

plot_chromatin_frac <- function(infile, outfile, width, height) {
    name <- gsub(".idxstats", "", basename(infile))
    p <- read.table(infile) %>%
        dplyr::select(1, 3) %>%
        rename_with(~ c("chrom", "counts")) %>%
        mutate(percent = round(counts * 100 / sum(counts), 2)) %>%
        dplyr::filter(chrom %in% paste0("chr", c(1:22, "X", "Y", "M"))) %>%
        mutate(chrom = factor(chrom, levels = rev(paste0("chr", c(1:22, "X", "Y", "M"))))) %>%
        ggplot(aes(percent, chrom)) +
        geom_col(aes(fill = chrom), show.legend = FALSE) +
        geom_text(aes(label = percent), hjust = -0.1) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
        labs(x = "Distribution of reads on different chromosomes(%)", y = "Chromosomes", title = name) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    pdf(outfile, width, height)
    print(p)
    dev.off()
}

# parse arguments
infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
width <- as.numeric(snakemake@params[["width"]])
height <- as.numeric(snakemake@params[["height"]])

cat("begin plot")
plot_chromatin_frac(infile, outfile, width, height)
cat("save to", outfile)
sink()