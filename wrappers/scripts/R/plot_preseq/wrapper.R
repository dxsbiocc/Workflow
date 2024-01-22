options(warn=-1)
# output to log file
sink(snakemake@log[[1]])

suppressWarnings(suppressMessages(library(tidyverse)))

plot_preseq <- function(infile, outfile, width, height) {
    name <- gsub(".lc_extrap", "", basename(infile))
    p <- suppressMessages(read_delim(infile, delim = '\t')) %>%
        rename_with(~c("total", "distinct", "lower", "upper")) %>%
        mutate_all(function(x) x / 1e6) %>% # per millions of reads
        ggplot(aes(x = total, y = distinct)) +
            geom_line(color = '#f1404b') +
            geom_ribbon(aes(ymin = lower, ymax = upper), fill = '#88dba3', alpha = 0.4) +
        labs(
            title = paste('Preseq estimated yield for', name),
            x = 'Sequenced fragments [ millions ]',
            y = 'Expected distinct fragments [ millions ]'
        ) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(outfile, p, width = width, height = height, device = "pdf")
}

# parse arguments
infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
width <- as.numeric(snakemake@params[["width"]])
height <- as.numeric(snakemake@params[["height"]])

cat("begin plot")
plot_preseq(infile, outfile, width, height)
cat("save to", outfile)
sink()