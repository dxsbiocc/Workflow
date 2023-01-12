# library
suppressMessages(library(ChIPseeker))
suppressMessages(library(GenomicFeatures))

annotate_peaks <- function(peakfile, gtf, sample, out_anno) {
  if (!dir.exists(out_anno))
      dir.create(out_anno)
  outfile <- file.path(out_anno, sample)

  suppressWarnings(tx <- makeTxDbFromGFF(file = gtf))
  peak <- readPeakFile(peakfile = peakfile, header = FALSE)
  suppressWarnings(peakAnno <- annotatePeak(peak = peak, TxDb = tx, assignGenomicAnnotation = TRUE))

  pdf(file = paste0(outfile, '.peakAnno.pdf'))
  plotAnnoPie(peakAnno, main = paste0(sample, "\nDistribution of Peaks"), line = -8)
  dev.off()

  write.table(as.data.frame(peakAnno@anno), file = paste0(outfile, '.peakAnno.txt'), sep = '\t', row.names = FALSE)
}

annotate_peaks(snakemake@input[['peakfile']], snakemake@input[['gtf']], snakemake@params[['sample']], snakemake@params[['output']])