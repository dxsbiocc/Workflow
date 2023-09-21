# Snakemake Workflow for NGS Analysis

## Introduction

The snakemake analysis workflow for bioinformatics analysis, including

- [x] ATAC-seq/Cut&Tag/ChIP-seq
- [X] DNA-seq(WGS and WES)
- [X] Assembly
- [ ] RNA-seq
- [ ] HiC

**To-Do list will be updated soon**

## Usage

#### **Configration**

1. One-Click Configuration

    ```sh
    sh setup.sh
    ```

2. Select the appropriate parameters in the file(`config.yaml`), such as reference genome.

3. Sample data information into tabular form.

#### **Sample data**

example files in directory `example`, you need to modify the path of the files `sample_info.json` and `sample_list.txt`

- *sample_list.txt*, use `generate_setting.py`

| sample | fastq1 | fastq2 | type(optional) |
| ------ | ------ | ------ | -------------- |
|   sp1  | path/to/sp1.R1.fq.gz | path/to/sp1.R2.fq.gz | chip |
|   sp2  | path/to/sp2.R1.fq.gz | path/to/sp2.R2.fq.gz | atac |

- *sample_info.json*

```json
{
    'sample1': 'control1',
    'sample2': 'control2'
}
```

#### **Running**

```sh
# run in local
snakemake -s path/to/Snakefile --use-conda -c4
# or run in the slurm task management system
snakemake -s path/to/Snakefile --profile path/to/config/slurm
```

## Description

### 1. ATAC-seq

- fastp: quality control and remove low quality reads
- bowtie2: mapping to reference genome
- samtools: sort and index bam file
- picard: remove duplicates
- deeptools/bedtools+samtools: filter reads with mapping quality < 30 or in blacklist region, mapped to mitochondria or unmapped reads.
    - deeptools: alignmentsieve, shift 9-bp [optional]
    - bedtools+samtools
- macs2: peak calling
- homer: annotation peaks and motif analysis
- deeptools: plot heatmap and profile

### 2. DNA-seq

- fastp: quality control and remove low quality reads
- bwa: mapping to reference genome
- samtools: sort and index bam file
- picard: remove duplicates
- samtools: statistics mapping information
- gatk: call variants
- annovar: annotate variants
- delly: detect SVs

### 3. Assembly

- hifiadapterfit: remove adapter for HIFI reads;
- fastp: qc for NGS reads(DNA-seq, RNA-seq, HiC);
- jellyfish: assessing heterozygosity of species, if DNA-seq exist;
- genomescope/genomescope2: plot k-mer histogram;
- hifiasm: sssembly(solo, hic, trio);
- gfatools/seqkit/quast: qc for assemble and get fasta;
- merqury/meryl: assessing heterozygosity for assemble.fasta;
- busco: sssessing gene assembly quality and completeness of gene predictions;
- juicer: HiC reads processing;
- 3d-dna: scaffolding using HiC reads(The /tmp directory requires a large space. If not, you need to specify `--tempdir` where parallel is used in the code);
- ...

## Notes

1. default running in conda env, most of the packages and software no need to install mannually
