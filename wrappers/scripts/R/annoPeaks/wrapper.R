sink(snakemake@log[[1]])
options(warn=-1)
# library
suppressMessages(library(ChIPseeker))
suppressMessages(library(GenomicFeatures))

annotate_peaks <- function(snakemake) {
  # parse parameters
  peakfile <- snakemake@input[["peakfile"]]
  gtf <- snakemake@input[["gtf"]]

  anno_txt <- snakemake@output[["anno_txt"]]
  anno_pdf <- snakemake@output[["anno_pdf"]]

  outdir <- dirname(anno_txt)
  if (!dir.exists(outdir))
    dir.create(outdir)

  if (file.size(peakfile) == 0) {
    # peak file empty
    file.create(anno_txt)
    file.create(anno_pdf)
    # print info to log file
    cat("The Peak file is empty!\n")
    return(1)
  }
  suppressWarnings(tx <- makeTxDbFromGFF(file = gtf))
  cat("Read peaks from file:", peakfile, "\n")
  peak <- suppressMessages(readPeakFile(peakfile = peakfile, header = FALSE))
  cat("Annotate peak to genes ...\n")
  peak_anno <- suppressMessages(annotatePeak(peak = peak, TxDb = tx, assignGenomicAnnotation = TRUE))

  cat("Plot the distribution of peaks.\n")
  pdf(file = anno_pdf)
  plotAnnoPie(peak_anno, main = paste0(snakemake@wildcards[["pair"]], "\nDistribution of Peaks"), line = -8)
  dev.off()

  cat("Save result to", anno_txt)
  write.table(as.data.frame(peak_anno@anno), file = anno_txt, sep = "\t", row.names = FALSE)
  cat("Annotation done!")
}

annotate_peaks(snakemake)
sink()