sink(snakemake@log[[1]])

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(ggrepel)))


plot_volcano <- function(data, pvalue = 0.05, lfc = 1, line.type = c("line", "curve"),
                         color = c("#2DB2EB", "#d8d8d8", "#EB4232"), facet = NULL, top = 5,
                         label.size = 10, top.label = NULL) {
    options(warn=-1)
    # lint type: line or curve
    line.type <- match.arg(line.type)
    if (line.type == "line") {
        data <- mutate(data, group = case_when(
            log2FC > lfc & p.value < pvalue ~ "Up",
            log2FC < -lfc & p.value < pvalue ~ "Down",
            TRUE ~ "None"
        ))
        my_line <- list(
            geom_vline(xintercept = c(-lfc, lfc), colour = "#737373", linetype = "dashed", linewidth = 0.6),
            geom_hline(yintercept = -log10(pvalue), colour = "#737373", linetype = "dashed", linewidth = 0.6)
        )
    } else if (line.type == "curve") {
        data %<>% mutate(
            curve = case_when(
                log2FC > 0 ~ 1 / (log2FC),
                log2FC <= 0 ~ 1 / (-log2FC)
            ),
            group = case_when(
                -log10(p.value) > curve & log2FC >= lfc ~ "Up",
                -log10(p.value) > curve & log2FC <= -lfc ~ "Down",
                TRUE ~ "None"
            )
        )
        my_fun <- function(left, start, right, pvalue, log2FC) {
            input1 <- c(seq(left, start, by = 0.01), start)
            input2 <- seq(start, right, by = 0.01)
            y1 <- 1 / abs(input1) + (-log10(pvalue))
            y2 <- 1 / (input2) + (-log10(pvalue))
            dff <- do.call(rbind, list(
                data.frame(x = input1 - log2FC, y = y1),
                data.frame(x = 0, y = NA),
                data.frame(x = input2 + log2FC, y = y2)
            ))

            return(dff)
        }
        my_line <- geom_line(
            data = my_fun(min(data$log2FC), min(data$p.value), max(data$log2FC), pvalue, lfc),
            mapping = aes(x, log10(y)), color = "black", linetype = "dashed", linewidth = 0.6
        )
    } else {
        stop("argument 'line.type' not support!")
    }
    if (is.null(top.label)) {
        if (!is.null(facet)) {
            top.label <- filter(data, group != "None") %>%
                group_by(!!sym(facet), group) %>%
                top_n(top, wt = abs(log2FC))
        } else {
            top.label <- filter(data, group != "None") %>%
                group_by(group) %>%
                top_n(top, abs(log2FC))
        }
    }
    # color
    mycolor <- structure(color, names = c("Down", "None", "Up"))
    color <- mycolor[sort(unique(data$group))]
    # plot
    p <- ggplot(data, aes(log2FC, -log10(p.value))) +
        geom_point(aes(colour = group), alpha = 0.7, size = 2.2) +
        my_line +
        geom_text_repel(
            data = top.label,
            mapping = aes(label = symbol),
            size = label.size, max.overlaps = 160,
            box.padding = unit(0.35, "lines"),
            point.padding = unit(0.3, "lines")
        ) +
        labs(title = " ", x = "log2(FoldChange)", y = "-log10(Pvalue)") +
        scale_colour_manual(values = color, name = NULL) +
        {
            if (!is.null(facet)) {
                facet_wrap(facets = facet, scales = "free_y")
            }
        } +
        theme_bw() +
        theme(
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 14),
            legend.text = element_text(size = 14),
            strip.background = element_rect(fill = "#b3de69"),
        )
    p
}

# parse arguments
infile <- snakemake@input[["all_genes"]]
outfile <- snakemake@output[["volcano"]]
logfile <- snakemake@log[[1]]

line_type <- snakemake@params[["line_type"]]
color <- unlist(strsplit(snakemake@params[["color"]], " "))

height <- as.numeric(snakemake@params[["height"]])
width <- as.numeric(snakemake@params[["width"]])

writeLines(paste("read plot data from", infile), logfile)
plot_data <- suppressMessages(read_csv(infile))
# main
writeLines(paste("volcano plot at", outfile), logfile)
pdf(outfile, width = width, height = height)
tryCatch(
    plot_volcano(plot_data, label.size = 3, top = 10, line.type = line_type, color = color),
    warning = function(w) {
        writeLines(paste("Warning: ", w$message), logfile)
  }
)
dev.off()
writeLines("plot volcano done!", logfile)

sink()