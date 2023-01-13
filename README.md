# Workflow

## Introduction

The snakemake analysis workflow for Epigenomic analysis, Such as

- ATAC-seq
- Cut&Tag
- Chip-seq

## Usage

#### Config

config file in `config/config.yaml`, Choose the right parameters.

#### Sample data

example files in `example`, you need to modify `sample_info.json` and `sample_list.txt`

#### Running

```sh
cd pipe
snakemake -s Snakefile --profile ../config/slurm
```

## Notes

1. change directory to `pipe`
2. default running in conda env, most of the packages and software no need to install mannually
3. python script `get_stats.py` need package `pandas`