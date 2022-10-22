# install packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("ChIPseeker", quietly = TRUE)) BiocManager::install("ChIPseeker")
if (!require("GenomicFeatures", quietly = TRUE)) BiocManager::install("GenomicFeatures")

# library
library(optparse)
library(ChIPseeker)
library(GenomicFeatures)

annotate_peaks <- function (peakfile, sample, gtf, out_anno) {
  if (!dir.exists(out_anno))
      dir.create(out_anno)
  outfile <- file.path(out_anno, sample)

  tx <- makeTxDbFromGFF(file = gtf)
  peak <- readPeakFile(peakfile = peakfile, header = FALSE)
  peakAnno <- annotatePeak(peak = peak, TxDb = tx, assignGenomicAnnotation = TRUE)

  pdf(file = paste0(outfile, '.peakNonotation.pdf'))
  plotAnnoPie(peakAnno, main = paste0(outfile, "\nDistribution of Peaks"), line = -8)
  dev.off()

  write.table(as.data.frame(peakAnno@anno), file = paste0(outfile, '.peakAnno.txt'), sep = '\t', row.names = FALSE)
}

annotate_peaks(nakemake@input[['peakfile']], nakemake@input[['sample']], nakemake@input[['gtf']], snakemake@output[[1]])