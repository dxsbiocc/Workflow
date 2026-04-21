# Materials and Methods

## 1. Overview of the Bioinformatics Analysis Framework

All next-generation sequencing (NGS) data analyses were performed using a comprehensive Snakemake-based workflow framework (version 7.12.0 or higher) developed in-house. This modular framework enables reproducible, scalable, and automated processing of ATAC-seq and RNA-seq data. The workflow leverages Conda environments for software dependency management and supports execution on both local workstations and high-performance computing clusters (SLURM).

## 2. ATAC-seq Data Processing and Analysis

### 2.1 Raw Data Preprocessing

Raw paired-end sequencing reads in FASTQ format were subjected to quality control and adapter trimming using **fastp** (version 0.23.2). The trimming process included:

- Removal of adapter sequences (Nextera or TruSeq adapters)
- Quality-based filtering with a minimum Phred quality score threshold
- Removal of reads with excessive N bases
- Automatic detection of read pairs

Alternative trimming tools available in the workflow include Trimmomatic, Cutadapt, and Trim Galore, selected based on specific experimental requirements.

### 2.2 Read Alignment

Trimmed reads were aligned to the reference genome (hg19/GRCh37 or GRCh38) using **Bowtie2** (version 2.5.1) with the following key parameters:

- `--very-sensitive` for high sensitivity alignment
- Proper handling of paired-end reads
- Maximum insert size considerations for ATAC-seq fragments

For experiments requiring different alignment strategies, alternative aligners including BWA-MEM, BWA-MEM2, HISAT2, and STAR were available within the workflow.

### 2.3 Post-Alignment Processing

#### 2.3.1 Duplicate Removal

PCR duplicates were removed using **Picard MarkDuplicates** (version 2.27.5) to eliminate technical duplicates arising from PCR amplification during library preparation. The duplicate removal process generated:

- De-duplicated BAM files for downstream analysis
- Duplicate metrics for quality assessment

Alternative duplicate removal was available via **Sambamba**.

#### 2.3.2 BAM Filtering and Quality Control

Post-duplication removal BAM files underwent stringent filtering:

1. **Alignment Sieve (deepTools)**: ATAC-seq specific coordinate shifting was applied using `alignmentSieve` with `--ATACshift` to account for the 9 bp offset of Tn5 transposase insertion.

2. **Blacklist Filtering**: Reads overlapping known blacklist regions (ENCODE blacklisted genomic regions) were removed using **BEDTools intersect**.

3. **Mitochondrial Read Removal**: Reads mapping to the mitochondrial chromosome (chrM) were excluded to reduce noise from mitochondrial DNA contamination.

4. **Mapping Quality Filtering**: Only reads with MAPQ ≥ 30 were retained for high-confidence analysis.

5. **Sorted and Indexed BAM Generation**: Filtered BAM files were coordinate-sorted and indexed using **SAMtools**.

### 2.4 Peak Calling

Open chromatin regions (peaks) were identified using **MACS2** (Model-based Analysis of ChIP-Seq, version 2.2.7.1) with experiment-specific parameters:

#### For ATAC-seq:
```
macs2 callpeak --bdg --SPMR --gsize hs --keep-dup all \
    --qvalue 0.05 --shift -75 --extsize 150 \
    --nomodel --call-summits --nolambda
```

#### For Cut&Tag:
```
macs2 callpeak --bdg --SPMR --gsize hs --keep-dup all \
    --qvalue 0.05 --shift -75 --extsize 150 \
    --nomodel --call-summits --nolambda
```

#### For ChIP-seq:
```
macs2 callpeak --bdg --SPMR --gsize hs --keep-dup all --qvalue 0.05
```

Two peak calling modes were supported:
- **Narrow peaks**: For transcription factor binding sites (default)
- **Broad peaks**: For histone modifications (e.g., H3K27me3, H3K36me3)

Peak calling was performed with matched input/IgG controls where available. The workflow supported paired-sample analysis (treatment vs. control) for accurate peak identification.

### 2.5 Peak Annotation and Quantification

#### 2.5.1 Peak Annotation

Called peaks were annotated to genomic features using the **ChIPseeker** R package, assigning peaks to:
- Promoter regions (±3 kb from transcription start sites)
- 5' UTR, 3' UTR
- Exons and introns
- Intergenic regions
- Distal intergenic regions

Annotation results included distribution plots and detailed text files mapping each peak to the nearest gene.

#### 2.5.2 Peak Quantification

Read counts within peak regions were quantified using **featureCounts** from the Subread package (version 2.0.6) with SAF (Simplified Annotation Format) annotation files derived from MACS2 peak calls. Parameters included:
- `-F SAF`: SAF format annotation
- Proper handling of paired-end reads
- Strand-specific counting when applicable

### 2.6 Motif Discovery

De novo motif discovery was performed using **HOMER** (Hypergeometric Optimization of Motif EnRichment, version 4.11) with the `findMotifsGenome.pl` command. Analysis parameters included:
- `-size 200`: 200 bp window around peak centers
- Motif enrichment against genomic background
- Known motif matching from HOMER database

### 2.7 Visualization and Signal Tracks

Normalized signal coverage tracks were generated using **deepTools bamCoverage** with:
- `--binSize 10`: 10 bp bin size
- `--normalizeUsing RPGC`: Reads Per Genomic Content normalization (1x coverage)
- Output in bigWig format for genome browser visualization

#### Profile and Heatmap Generation

deepTools was used to generate:
- **TSS profiles**: Signal distribution around transcription start sites (±3 kb)
- **Gene body profiles**: Signal across gene bodies with 5 kb extension
- **Heatmaps**: Clustered heatmaps of signal intensity
- **Correlation analysis**: Spearman correlation of read counts between samples

### 2.8 Quality Control Metrics

Comprehensive QC was performed following ENCODE ATAC-seq pipeline standards:

#### 2.8.1 Library Complexity Metrics
- **PBC1 (PCR Bottleneck Coefficient 1)**: M1/MDISTINCT
  - M1: Number of genomic locations with exactly one unique read
  - MDISTINCT: Number of distinct genomic locations with unique reads
- **PBC2 (PCR Bottlenecking Coefficient 2)**: M1/M2
- **NRF (Non-Redundant Fraction)**: Distinct uniquely mapping reads / Total reads

Quality thresholds (ATAC-seq):
- Ideal: PBC1 ≥ 0.9, PBC2 > 3, NRF > 0.9
- Acceptable: 0.7 ≤ PBC1 ≤ 0.9, 1 ≤ PBC2 ≤ 3, 0.7 ≤ NRF ≤ 0.9
- Concerning: PBC1 < 0.7, PBC2 < 1, NRF < 0.7

#### 2.8.2 Enrichment Metrics
- **FRiP (Fraction of Reads in Peaks)**: Proportion of reads falling within called peaks (ideal > 0.3)
- **Fraction of reads in annotated regions**: Distribution across DHS, promoters, enhancers
- **Chromosome coverage distribution**: Read count distribution across chromosomes

#### 2.8.3 Insert Size Distribution
Fragment length distribution was analyzed using **Picard CollectInsertSizeMetrics** with histogram width of 1000 bp and maximum of 5 million reads for efficiency. ATAC-seq fragment sizes were expected to show:
- Nucleosome-free regions (<100 bp)
- Mono-nucleosome (~200 bp)
- Di-nucleosome (~400 bp)

#### 2.8.4 GC Bias Analysis
GC bias metrics were collected using **Picard CollectGcBiasMetrics** to assess library preparation artifacts.

#### 2.8.5 Library Complexity Prediction
**Preseq lc_extrap** was used to predict library complexity and estimate the depth of unique molecules at higher sequencing depths.

#### 2.8.6 Summary Report
All QC metrics were aggregated into a comprehensive Excel summary report for quality assessment and documentation.

### 2.9 Spike-in Normalization (Optional)

For experiments requiring absolute quantification, *E. coli* spike-in normalization was available. Reads were mapped to the *E. coli* genome separately, and normalization factors were calculated based on spike-in read proportions.

### 2.10 Down-sampling (Optional)

For comparative analyses, samples were down-sampled to equal read depths using **SAMtools view** with random subsampling to ensure fair comparison between samples.

## 3. RNA-seq Data Processing and Analysis

### 3.1 Raw Data Preprocessing

RNA-seq raw reads were processed similarly to ATAC-seq using **fastp** for adapter trimming and quality filtering. Key considerations included:

- Stranded library preservation
- Quality score thresholding
- Adapter sequence removal (Illumina TruSeq adapters)

### 3.2 Read Alignment

#### 3.2.1 Splice-Aware Alignment

For gene-level and transcript-level quantification requiring alignment:

**STAR** (Spliced Transcripts Alignment to a Reference, version 2.7.10b) was used with:
- `--sjdbOverhang ReadLength-1`: Optimized for 150 bp reads
- `--sjdbGTFfile`: Gene annotation-guided spliced alignment
- `--quantTranscriptomeBan IndelSoftclipSingleend`: RSEM-compatible output
- `--outSAMattrRGline`: Read group information for sample tracking

**HISAT2** (Hierarchical Indexing for Spliced Alignment of Transcripts 2) was available as an alternative splice-aware aligner.

#### 3.2.2 Splice-Unaware Alignment

For certain quantification methods (featureCounts, HTSeq), **Bowtie2** or **BWA-MEM** were used.

### 3.3 Post-Alignment Processing

Aligned reads were processed through duplicate removal using **Picard MarkDuplicates** (similar to ATAC-seq pipeline), followed by coordinate sorting and indexing using **SAMtools**.

### 3.4 Gene and Transcript Quantification

The workflow supported multiple quantification strategies:

#### 3.4.1 RSEM (RNA-Seq by Expectation-Maximization)

**RSEM** (version 1.3.3) was used for accurate transcript-level quantification:

1. **Reference Preparation**: RSEM reference indices were prepared from the genome FASTA and GTF annotation using `rsem-prepare-reference`.

2. **Quantification Modes**:
   - FASTQ mode: Direct quantification from trimmed reads with embedded alignment
   - BAM mode: Quantification from pre-aligned BAM files

3. **Parameters**:
   - `--estimate-rspd`: Estimate read start position distribution
   - `--calc-ci`: Calculate 95% credibility intervals
   - `--seed 123456`: Random seed for reproducibility
   - `--no-bam-output`: Suppress BAM output for space efficiency

4. **Outputs**: Gene-level and isoform-level counts, TPM, and FPKM values

#### 3.4.2 Salmon (Selective Alignment)

**Salmon** (version 1.10.0) provided rapid transcript-level quantification using selective alignment:

1. **Index Building**: Transcriptome index was built using `salmon index`

2. **Quantification**:
   - Library type auto-detection (`-l A`)
   - Bias correction algorithms enabled
   - Output in quant.gene.sf format

3. **Advantages**: Pseudo-alignment approach for fast, accurate quantification without full alignment

#### 3.4.3 Kallisto (Pseudo-alignment)

**Kallisto** (version 0.48.0) was available for ultra-fast transcript quantification via pseudo-alignment using k-mer matching.

#### 3.4.4 featureCounts (Read Counting)

**featureCounts** from the Subread package provided gene-level read counting from aligned BAM files:

```
featureCounts -O --fracOverlap 0.2 -J -g gene_id -a annotation.gtf
```

Parameters:
- `-O`: Allow reads to overlap multiple features
- `--fracOverlap 0.2`: Minimum 20% fractional overlap
- `-J`: Count junction reads
- `-g gene_id`: Group by gene_id attribute
- Strand-specific counting when applicable

#### 3.4.5 HTSeq-count

**HTSeq** was available as an alternative read counting tool with intersection-nonempty mode for handling overlapping features.

### 3.5 Expression Matrix Generation

Sample-level quantification results were merged into a unified expression matrix using custom Python scripts. The matrix included:
- Raw read counts for differential expression analysis
- TPM (Transcripts Per Million) normalized values
- FPKM (Fragments Per Kilobase Million) values when applicable

**Normalization Formulas**:

TPM (Transcripts Per Million):
$$TPM_i = \frac{\frac{q_i}{l_i}}{\sum_j(\frac{q_j}{l_j})} \times 10^6$$

where $q_i$ is the read count for transcript $i$, and $l_i$ is the effective length.

FPKM (Fragments Per Kilobase Million):
$$FPKM_i = \frac{q_i}{l_i \times \sum_j q_j} \times 10^9$$

### 3.6 Differential Expression Analysis

Differential gene expression analysis was performed using multiple R/Bioconductor packages:

#### 3.6.1 DESeq2 (Default)

**DESeq2** (version 1.38.0) was the primary tool for differential expression:

1. **Input**: Raw count matrix and sample metadata
2. **Normalization**: Size factor normalization for library depth differences
3. **Dispersion Estimation**: Empirical Bayes shrinkage of dispersion estimates
4. **Statistical Testing**: Wald test with Benjamini-Hochberg FDR correction
5. **Filtering Criteria**:
   - Adjusted p-value < 0.05
   - |log2FoldChange| ≥ 1

#### 3.6.2 edgeR

**edgeR** (Empirical Analysis of Digital Gene Expression Data in R) was available using:
- TMM (Trimmed Mean of M-values) normalization
- Negative binomial generalized linear models
- Empirical Bayes quasi-likelihood F-tests

#### 3.6.3 limma-voom

**limma** (Linear Models for Microarray Data) with voom transformation was available for RNA-seq data, particularly useful for small sample sizes.

### 3.7 Visualization of Differential Expression Results

#### 3.7.1 Volcano Plots

Volcano plots were generated using ggplot2 to visualize:
- Statistical significance (-log10 adjusted p-value) on y-axis
- Magnitude of change (log2 fold change) on x-axis
- Color-coded points: upregulated (red), downregulated (blue), non-significant (gray)

#### 3.7.2 Principal Component Analysis (PCA)

PCA was performed on variance-stabilized transformed counts to:
- Assess sample clustering and batch effects
- Visualize relationships between biological replicates
- Identify outliers

### 3.8 Gene Fusion Detection (Optional)

For cancer-related studies, gene fusion detection was supported:

#### 3.8.1 Arriba

**Arriba** (version 2.4.0) detected gene fusions from STAR-aligned BAM files:
- Utilized STAR chimeric alignments
- Blacklist filtering for known artifacts
- Visualization of fusion breakpoints
- Optional: Structural variant file integration

#### 3.8.2 STAR-Fusion

**STAR-Fusion** (version 1.12.0) provided alternative fusion detection:
- CTAT genome library for annotation
- Kick-start mode from STAR chimeric junctions
- Standalone mode from FASTQ files

### 3.9 Functional Enrichment Analysis

Differentially expressed genes were subjected to pathway and functional enrichment analysis using **clusterProfiler**:

#### 3.9.1 Gene Ontology (GO) Enrichment

- Biological Process (BP)
- Molecular Function (MF)
- Cellular Component (CC)

#### 3.9.2 Pathway Enrichment

Supported pathway databases:
- **KEGG** (Kyoto Encyclopedia of Genes and Genomes)
- **MSigDB** (Molecular Signatures Database)
- **Reactome**

Parameters:
- p-value cutoff: 0.05
- q-value cutoff: 0.05 (FDR corrected)
- Species-specific annotation databases (org.Hs.eg.db for human)

## 4. Common Tools and Resources

### 4.1 Reference Genomes and Annotations

- **Human Genome**: GRCh37/hg19 or GRCh38/GENCODE v42
- **Gene Annotation**: GENCODE or RefGene GTF files
- **Blacklist Regions**: ENCODE blacklisted genomic regions

### 4.2 Software Versions and Environment Management

All software dependencies were managed through Conda environments defined in the workflow. Key tools included:

| Tool | Version | Purpose |
|------|---------|---------|
| Snakemake | ≥7.12.0 | Workflow orchestration |
| fastp | 0.23.2 | Read trimming |
| Bowtie2 | 2.5.1 | Read alignment |
| STAR | 2.7.10b | Splice-aware alignment |
| BWA-MEM | 0.7.17 | Alignment |
| Picard | 2.27.5 | Duplicate removal, QC |
| SAMtools | 1.17 | BAM manipulation |
| BEDTools | 2.30.0 | Genome arithmetic |
| MACS2 | 2.2.7.1 | Peak calling |
| deepTools | 3.5.1 | Visualization, QC |
| HOMER | 4.11 | Motif analysis |
| RSEM | 1.3.3 | Transcript quantification |
| Salmon | 1.10.0 | Pseudo-alignment |
| Kallisto | 0.48.0 | Pseudo-alignment |
| featureCounts | 2.0.6 | Read counting |
| HTSeq | 2.0.2 | Read counting |
| DESeq2 | 1.38.0 | Differential expression |
| edgeR | 3.40.0 | Differential expression |
| limma | 3.54.0 | Differential expression |
| clusterProfiler | 4.6.0 | Enrichment analysis |
| Arriba | 2.4.0 | Fusion detection |
| STAR-Fusion | 1.12.0 | Fusion detection |
| Preseq | 2.0.3 | Library complexity |
| MultiQC | 1.14 | Report aggregation |

### 4.3 Workflow Execution

The workflow was executed using:
```bash
snakemake -s path/to/Snakefile --use-conda -c 32
```

For cluster execution (SLURM):
```bash
snakemake -s path/to/Snakefile --profile path/to/config/slurm
```

### 4.4 Quality Control Summary

All analyses included comprehensive quality control at each step:
- Raw read quality (FastQC)
- Trimming statistics
- Alignment rates and uniqueness
- Duplicate rates
- Library complexity
- Signal-to-noise ratios

QC results were aggregated using **MultiQC** for comprehensive reporting.

## 5. Data Availability and Reproducibility

All raw sequencing data have been deposited in [appropriate repository, e.g., NCBI GEO/SRA, ENA]. The complete Snakemake workflow code, configuration files, and processed data are available upon reasonable request. The workflow ensures full reproducibility through:
- Version-controlled software environments
- Conda environment specifications
- Random seed setting for stochastic processes
- Comprehensive logging of all analysis steps

---

*Note: This Materials and Methods section describes the comprehensive bioinformatics framework implemented in this study. Specific parameter choices and tool selections may vary based on experimental design and research questions. All analyses were performed following established best practices in the field and ENCODE consortium guidelines where applicable.*
