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

def get_paired_bam(wildcards):
    paired = {'treatment': opj(OUTDIR, 'dedup/{pair}/{pair}.filtered.bam')}
    if SAMPLE_MAP and SAMPLE_MAP.get(wildcards.pair):
        paired['control'] = opj(OUTDIR, "dedup/{sp}/{sp}.filtered.bam").format(sp=SAMPLE_MAP[wildcards.pair])
    return paired
    
def get_bigwig(wildcards):
    sp_list = [opj(OUTDIR, "macs2/bigwig/{}.norm.bw").format(wildcards.control)]
    if CONTROL.get(wildcards.control):
        sp_list.extend([opj(OUTDIR, "macs2/bigwig/{}.norm.bw").format(sp) for sp in CONTROL[wildcards.control]])
    return sp_list

def get_macs2(sample, narrow=True):
    if not narrow:
        return MACS2_BAM + " -g hs --broad --broad-cutoff 0.05 -B --SPMR --keep-dup all"
    if MACS2_MAP:
        if 'atac' in MACS2_MAP[sample].lower():
            return MACS2_BAM + config['parameters']['macs2']['atac']
        elif 'chip' in MACS2_MAP[sample].lower():
            return MACS2_BAM + config['parameters']['macs2']['chip']
    else:
        return MACS2_BAM + config['parameters']['macs2']['atac']