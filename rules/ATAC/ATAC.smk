from collections import defaultdict


# macs2 parameters
MACS2_MAP = {}
if 'type' in DATA.columns:
    MACS2_MAP = DATA['type'].to_dict()
# downsample to the same depth
DOWNSAMPLE = config['control'].get('downsample')
if DOWNSAMPLE and isinstance(DOWNSAMPLE, int):
    OUTPUT.append(expand(opj(OUTDIR, "dedup/{sample}/{sample}.sampled.bam"), sample=SAMPLES))
# sample informations
if SAMPLE_MAP:
    CONTROL = defaultdict(list)
    for k, v in SAMPLE_MAP.items():
        if v:
            CONTROL[v].append(k)
        else:
            CONTROL[k] = []
else:
    CONTROL = {k : '' for k in PAIRS}

# macs
if config['control']['paired']:
    MACS2_BAM = '-f BAMPE '
else:
    MACS2_BAM = '-f BAM '
# peak mode
PEAKMODE = config['control']['peak_mode'].lower()
assert PEAKMODE in ['narrow', 'broad'], "`peak_mode` must in ['narrow', 'broad']"
if PEAKMODE == 'broad':
    OUTPUT.append(expand(opj(OUTDIR, "macs2/broad/{pair}_peaks.broadPeak"), pair=PAIRS))
else:
    OUTPUT.append(expand(opj(OUTDIR, "macs2/narrow/{pair}_peaks.narrowPeak"), pair=PAIRS))
############################################################
#                         0. Include                       #
############################################################
include: os.path.join(PATH, "rules/ATAC/common.smk")
include: os.path.join(PATH, "rules/common/trimmed.smk")
include: os.path.join(PATH, "rules/common/mapping.smk")
include: os.path.join(PATH, "rules/common/dedup.smk")
include: os.path.join(PATH, "rules/common/stats.smk")
include: os.path.join(PATH, "rules/ATAC/spikein.smk")
include: os.path.join(PATH, "rules/ATAC/filtering.smk")
include: os.path.join(PATH, "rules/ATAC/callpeak.smk")
include: os.path.join(PATH, "rules/ATAC/motif.smk")
include: os.path.join(PATH, "rules/ATAC/quantify.smk")
include: os.path.join(PATH, "rules/ATAC/quality_control.smk")
include: os.path.join(PATH, "rules/ATAC/report.smk")

############################################################
#                           Runing                         #
############################################################
# report: "report/workflow.rst"
# local
localrules:
    mapped_stats, mapped_plotBamStats, mapped_idxstats, mapped_flagstat, \
    dedup_stats, dedup_plotBamStats, dedup_idxstats, dedup_flagstat, \
    rmdup_index, get_peak, plot_chromosome_coverage, plot_preseq, plot_insertsizemetrics, \
    plot_gcbiasmetrics, frip, preseq, frac_anno_regions, collectinsertsizemetrics, qc_summary

rule use_all:
    input:
        # trimmed adapter
        expand(opj(OUTDIR, "trimmed/{sample}/{sample}.clean.{unit}.fq.gz"), sample=SAMPLES, unit=RUN),
        # spikein
        expand(opj(OUTDIR, "mapped/{sample}/{sample}.spikein.bam"), sample=SAMPLES),
        # mapping and remove duplicates
        expand(opj(OUTDIR, "dedup/{sample}/{sample}.rmdup.bam.bai"), sample=SAMPLES),
        # shift 9bp
        expand(opj(OUTDIR, "dedup/{sample}/{sample}.filtered.bam"), sample=SAMPLES),
        # convert to bigwig
        expand(opj(OUTDIR, "macs2/bigwig/{sample}.norm.bw"), sample=SAMPLES),
        # motify analysis
        expand(opj(OUTDIR, "motifs/{pair}"), pair=PAIRS),
        # peak annotate
        expand(opj(OUTDIR, 'macs2/anno/{pair}.peakAnno.{ext}'), pair=PAIRS, ext=['pdf', 'txt']),
        # bw data quality
        expand(opj(OUTDIR, "macs2/matrix/{control}.{pos}.{types}.pdf"), control=CONTROL.keys(), pos=['tss', 'genebody'], types=['heatmap', 'profile']),
        # correlation heatmap
        opj(OUTDIR, "macs2/bigwig/heatmap_spearman_corr_readCounts.pdf"),
        # QC
        # stats
        expand(opj(OUTDIR, "QC/{bam}/stats/{sample}.{stats}"), bam=["mapped", "dedup"], sample=SAMPLES, stats=['stats', 'idxstats', 'flagstats']),
        expand(opj(OUTDIR, "QC/{bam}/plot/{sample}"), bam=["mapped", "dedup"], sample=SAMPLES),
        expand(opj(OUTDIR, "QC/chromosome_coverage/{sample}.pdf"), sample=SAMPLES),
        expand(opj(OUTDIR, "QC/frip/{pair}.txt"), pair=PAIRS),
        expand(opj(OUTDIR, "QC/preseq/{sample}.pdf"), sample=SAMPLES),
        expand(opj(OUTDIR, "QC/anno_regions/{sample}.txt"), sample=SAMPLES),
        expand(opj(OUTDIR, "QC/insertsize/{sample}.isize_dist.pdf"), sample=SAMPLES),
        expand(opj(OUTDIR, "QC/gcbias/{sample}.fraglen_dist.pdf"), sample=SAMPLES),
        expand(opj(OUTDIR, "QC/lib_complexity/{sample}.txt"), sample=SAMPLES),
        opj(OUTDIR, "QC/QC_Summary.xlsx")
    message:
        print('ATAC-seq pipeline is runing!')