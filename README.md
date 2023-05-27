# Workflow

## Introduction

The snakemake analysis workflow for bioinformatics analysis, including

- [x] ATAC-seq/Cut&Tag/ChIP-seq
- [ ] RNA-seq
- [ ] DNA-seq(WGS and WES)
- [ ] HiC

**To-Do list will be updated soon**

## Usage

#### **Configration**

1. install packages

    ```sh
    pip install -r requirement.txt
    ```

2. copy the file `base.py` in `utils` to package `snakemake-wrapper-utils` installed directory.
3. setting the config file `config.yaml` in `config`, choose the properly parameters.
4. setting the `root_dir` value in `config/config.yaml`

#### **Sample data**

example files in directory `example`, you need to modify the path of the files `sample_info.json` and `sample_list.txt`

*sample_list.txt*

| sample | fastq1 | fastq2 | type(optional) |
| ------ | ------ | ------ | -------------- |
|   sp1  | path/to/sp1.R1.fq.gz | path/to/sp1.R2.fq.gz | chip |
|   sp2  | path/to/sp2.R1.fq.gz | path/to/sp2.R2.fq.gz | atac |

*sample_info.json*

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
