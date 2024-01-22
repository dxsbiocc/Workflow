options(warn = -1)
suppressWarnings(suppressMessages(library(tidyverse)))

plot_gc <- function(infile, outfile, width, height) {
    name <- gsub(".gcmetrics.txt", "", basename(infile))
    p <- suppressMessages(read_delim(infile, comment = "#")) %>% ggplot(aes(x = GC)) +
        geom_line(aes(y = NORMALIZED_COVERAGE, color = "Normalized coverage")) +
        geom_line(aes(y = MEAN_BASE_QUALITY, color = "Mean base quality at GC%")) +
        geom_line(aes(y = WINDOWS / sum(WINDOWS), color = "Windows at GC%")) +
        scale_x_continuous(limits = c(0, 100)) +
        labs(y = "Normalized coverage", color = "Legend", title = name) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5))

    ggsave(outfile, p, device = "pdf", width = width, height = height)
}

# parse arguments
infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
width <- as.numeric(snakemake@params[["width"]])
height <- as.numeric(snakemake@params[["height"]])

plot_gc(infile, outfile, width, height)
