sink(snakemake@log[[1]])
# library
suppressMessages(library(ChIPseeker))
suppressMessages(library(GenomicFeatures))

annotate_peaks <- function(snakemake) {
  # parse parameters
  peakfile <- snakemake@input[["peakfile"]]
  gtf <- snakemake@input[["gtf"]]
  sample <- snakemake@params[["sample"]]
  out_anno <- snakemake@params[["output"]]
  if (!dir.exists(out_anno))
    dir.create(out_anno)
  outfile <- file.path(out_anno, sample)

  if (file.size(peakfile) == 0) {
    # peak file empty
    file.create(paste0(outfile, ".peakAnno.txt"))
    file.create(paste0(outfile, ".peakAnno.pdf"))
    # print info to log file
    writeLines("The Peak file is empty!", snakemake@log)
    return(1)
  }
  suppressWarnings(tx <- makeTxDbFromGFF(file = gtf))
  peak <- readPeakFile(peakfile = peakfile, header = FALSE)
  peak_anno <- annotatePeak(peak = peak, TxDb = tx, assignGenomicAnnotation = TRUE)

  pdf(file = paste0(outfile, ".peakAnno.pdf"))
  plotAnnoPie(peak_anno, main = paste0(sample, "\nDistribution of Peaks"), line = -8)
  dev.off()

  write.table(as.data.frame(peak_anno@anno), file = paste0(outfile, ".peakAnno.txt"), sep = "\t", row.names = FALSE)

  writeLines("Annotation done!", snakemake@log)
}

annotate_peaks(snakemake)
sink()