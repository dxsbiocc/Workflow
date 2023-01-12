# Workflow

## introduction

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
snakemake -s Snakefile --profile ../config/slurm
```