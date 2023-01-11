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
# database path
PATH = os.path.abspath(os.getcwd())
if config['data']['db'] == 'local':
    DATABASE = os.path.join(PATH, '../data')
else:
    DATABASE = config['data']['db']
# blacklist
BLACKLIST = os.path.join(DATABASE, 'blacklist', f'{genome}.blacklist.bed')
# GC
GC = os.path.join(DATABASE, 'GC', f'{genome}.gc')
# paired
if config['control']['paired']:
    RUN = ['R1', 'R2']
else:
    RUN = ['R1']

SPIKEIN = True if config['control']['spike_in'] else False
SHIFT = True if config['control']['shift'] else False
COLORS = ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6"]

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

def get_control_bam(wildcards):
    if SAMPLE_MAP:
        sp = "dedup/{}.filtered.bam".format(SAMPLE_MAP[wildcards.pair])
    else:
        sp = ""
    return sp
    
def get_bigwig(wildcards):
    if SAMPLE_MAP:
        sp_list = ["macs2/bigwig/{}.norm.bw".format(wildcards.pair), "macs2/bigwig/{}.norm.bw".format(SAMPLE_MAP[wildcards.pair])]
    else:
        sp_list = "macs2/bigwig/{}.norm.bw".format(wildcards.pair)
    return sp_list

def get_wrapper(*args, local=True):
    """Get wrappers path"""
    if local:
        abspath = os.path.join(PATH, '../wrappers', *args)
        return "file://{}".format(abspath)
    else:
        raise ValueError("Please use local version!")

def get_script(script):
    return os.path.join(PATH, '../scripts', script)

def get_adapter(method='fastp'):
    """get adapter path"""
    if config['control']['adapters'] == 'Truseq':
        adapt_file = os.path.join(DATABASE, "adapters", "Truseq3.PE.fa")
    elif config['control']['adapters'] == 'Nextera':
        adapt_file = os.path.join(DATABASE, "adapters", "NexteraPE-PE.fa")
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
