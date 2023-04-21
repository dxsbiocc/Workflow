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
BLACKLIST = os.path.join(DATABASE, 'blacklist', f'{genome}.blacklist.bed')
# GC
GC = os.path.join(DATABASE, 'GC', f'{genome}.gc')
COLORS = ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6"]

# control parameters
SPIKEIN = True if config['control']['spike_in'] else False
SHIFT = True if config['control']['shift'] else False
# -------------------- Helper functions -------------------- #
# def get_fastq(wildcards):
#     """Get fastq files of given sample-unit."""
#     fastqs = samples.loc[wildcards.sample, ["fastq1", "fastq2"]].dropna()
#     if len(fastqs) == 2:
#         return {'fq1': fastqs.fastq1, 'fq2': fastqs.fastq2}
#     return {'fq1': fastqs.fastq1}
def get_control_bam(wildcards):
    if SAMPLE_MAP:
        sp = "dedup/{}.filtered.bam".format(SAMPLE_MAP[wildcards.pair])
    else:
        sp = ""
    return sp
    
def get_bigwig(wildcards):
    if SAMPLE_MAP:
        sp_list = [
            "macs2/bigwig/{}.norm.bw".format(wildcards.control), 
            *["macs2/bigwig/{}.norm.bw".format(sp) for sp in CONTROL[wildcards.control]]
        ]
    else:
        sp_list = "macs2/bigwig/{}.norm.bw".format(wildcards.control)
    return sp_list