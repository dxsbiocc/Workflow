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
# chromosome
CHROM_BED = os.path.join(DATABASE, 'chromsize', f'{genome}.chromsize.bed')
# blacklist
BLACKLIST = os.path.join(DATABASE, 'blacklist', f'{genome}.blacklist.bed')
# GC
GC = os.path.join(DATABASE, 'GC', f'{genome}.gc')
COLORS = ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6"]

# control parameters
SHIFT = True if config['control']['shift'] else False
# -------------------- Helper functions -------------------- #
# used for callpeak
def get_paired_bam(wildcards):
    suffix = '.filtered.bam'
    if DOWNSAMPLE:
        suffix = f'.sampled.bam'
    paired = {'treatment': opj(OUTDIR, f"dedup/{wildcards.pair}/{wildcards.pair}{sufix}")}
    if SAMPLE_MAP and SAMPLE_MAP.get(wildcards.pair):
        paired['control'] = opj(OUTDIR, "dedup/{sp}/{sp}" + suffix).format(sp=SAMPLE_MAP[wildcards.pair])
    return paired
# used for computematrix
def get_bigwig(wildcards):
    sp_list = [opj(OUTDIR, "macs2/bigwig/{}.norm.bw").format(wildcards.control)]
    if CONTROL.get(wildcards.control):
        sp_list.extend([opj(OUTDIR, "macs2/bigwig/{}.norm.bw").format(sp) for sp in CONTROL[wildcards.control]])
    return sp_list
# get macs2 parameters
def get_macs2(sample, narrow=True):
    if not narrow:
        return MACS2_BAM + " -g hs --broad --broad-cutoff 0.05 -B --SPMR --keep-dup all"
    if MACS2_MAP:
        if 'atac' in MACS2_MAP[sample].lower():
            return MACS2_BAM + config['parameters']['macs2']['callpeak']['atac']
        elif 'cut&tag' in MACS2_MAP[sample].lower():
            return MACS2_BAM + config['parameters']['macs2']['callpeak']['cutntag']
        elif 'chip' in MACS2_MAP[sample].lower():
            return MACS2_BAM + config['parameters']['macs2']['callpeak']['chip']
    else:
        return MACS2_BAM + config['parameters']['macs2']['callpeak']['cutntag']

# used for FRIP
def get_peakfile(wildcards):
    if PEAKMODE == 'broad':
        peaks = opj(OUTDIR, f"macs2/broad/{wildcards.pair}_peaks.broadPeak")
    else:
        peaks = opj(OUTDIR, f"macs2/narrow/{wildcards.pair}_peaks.narrowPeak")
    return peaks