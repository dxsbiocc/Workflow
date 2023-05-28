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

2. Choose pipeline and set the file path in `config.yaml`, such as reference genome.

3. sample information to table or json.

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
    'sp1': 'control_sample',
    'sp2': 'control_sample'
}
```

#### **Running**

```sh
cd pipe
# run in local
snakemake -s Snakefile --use-conda -c4
# or run in the slurm task management system
snakemake -s Snakefile --profile ../config/slurm
```

## Notes

1. change directory to `pipe`
2. default running in conda env, most of the packages and software no need to install mannually
3. python script `get_stats.py` need package `pandas`
