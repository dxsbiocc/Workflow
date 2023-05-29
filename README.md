# Workflow

## Introduction

The snakemake analysis workflow for bioinformatics analysis, including

- [x] ATAC-seq/Cut&Tag/ChIP-seq
- [X] DNA-seq(WGS and WES)
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

## Notes

1. change directory to `pipe`
2. default running in conda env, most of the packages and software no need to install mannually
3. python script `get_stats.py` need package `pandas`
