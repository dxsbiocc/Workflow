import os
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version


# -------------------- Other Parameters -------------------- #
genome_size = {
    'hg19': 2864785220,
    'hg38': 2913022398,
    'mm10': 2652783500,
    'mm9': 2620345972,
}
genome = config['data']['genome']
# chrom size
CHROM_SIZE = genome_size[genome]
# blacklist
BLACKLIST = os.path.join(config['data']['db'], 'blacklist', f'{genome}.blacklist.bed')
# GC
GC = os.path.join(config['data']['db'], 'GC', f'{genome}.gc')

# -------------------- Helper functions -------------------- #
# def get_fastq(wildcards):
#     """Get fastq files of given sample-unit."""
#     fastqs = samples.loc[wildcards.sample, ["fastq1", "fastq2"]].dropna()
#     if len(fastqs) == 2:
#         return {'fq1': fastqs.fastq1, 'fq2': fastqs.fastq2}
#     return {'fq1': fastqs.fastq1}

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = samples.loc[wildcards.sample, ["fastq1", "fastq2"]].dropna()
    if len(fastqs) == 2:
        return [fastqs.fastq1, fastqs.fastq2]
    return [fastqs.fastq1]

def get_wrapper(*args):
    """Get wrappers path"""
    return os.path.join(config['wrappers'], *args)

def get_script(script):
    return os.path.join(config['scripts'], script)

def get_adapter(method='fastp'):
    """get adapter path"""
    if config['control']['adapters'] == 'Truseq':
        adapt_file = os.path.join(config["data"]["db"], "adapters", "Truseq3.PE.fa")
    elif config['control']['adapters'] == 'Nextera':
        adapt_file = os.path.join(config["data"]["db"], "adapters", "NexteraPE-PE.fa")
    else:
        raise ValueError('%s not support!' % config['control']['adapters'])

    if method == 'fastp':
        pat = f"--adapter_fasta {adapt_file}"
    elif method == 'trimmomatic':
        pat = f'ILLUMINACLIP:{adapt_file}:2:15:4:4:true'
    elif method == 'cutadapt':
        pat = f'-a "file:{adapt_file}"'
    else:
        raise ValueError(f'{method} not support!')
    return pat
