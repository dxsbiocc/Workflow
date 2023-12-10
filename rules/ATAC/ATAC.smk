from collections import defaultdict

# configfile
configfile: os.path.join(PATH, "config/ATAC-seq.yaml")

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

if config['control']['paired']:
    MACS2_BAM = '-f BAMPE '
else:
    MACS2_BAM = '-f BAM '

# macs2 parameters
MACS2_MAP = {}
if 'type' in DATA.columns:
    MACS2_MAP = DATA['type'].to_dict()
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
include: os.path.join(PATH, "rules/ATAC/quality.smk")
include: os.path.join(PATH, "rules/ATAC/report.smk")

############################################################
#                           Runing                         #
############################################################
# report: "report/workflow.rst"
# local
localrules:
    stats, idxstats, flagstat

rule use_all:
    input:
        # trimmed adapter
        expand(opj(OUTDIR, "trimmed/{sample}/{sample}.clean.{unit}.fq.gz"), sample=SAMPLES, unit=RUN),
        # spikein
        expand(opj(OUTDIR, "mapped/{sample}/{sample}.spikein.bam"), sample=SAMPLES),
        # stats
        expand(opj(OUTDIR, "report/{bam}/stats/{sample}.{stats}"), bam=["mapped", "dedup"], sample=SAMPLES, stats=['stats', 'idxstats', 'flagstats']),
        expand(opj(OUTDIR, "report/{bam}/plot/{sample}"), bam=["mapped", "dedup"], sample=SAMPLES),
        # mapping and remove duplicates
        expand(opj(OUTDIR, "dedup/{sample}/{sample}.rmdup.bam.bai"), sample=SAMPLES),
        # shift 9bp
        expand(opj(OUTDIR, "dedup/{sample}/{sample}.filtered.bam"), sample=SAMPLES),
        # callpeak
        expand(opj(OUTDIR, "macs2/narrow/{pair}_peaks.narrowPeak"), pair=PAIRS),
        expand(opj(OUTDIR, "macs2/broad/{pair}_peaks.broadPeak"), pair=PAIRS),
        # convert to bigwig
        expand(opj(OUTDIR, "macs2/bigwig/{sample}.norm.bw"), sample=SAMPLES),
        # extract peaks
        expand(opj(OUTDIR, "macs2/narrow/{pair}_homer_peaks.txt"), pair=PAIRS),
        # motify analysis
        expand(opj(OUTDIR, "motifs/{pair}"), pair=PAIRS),
        # peak annotate
        expand(opj(OUTDIR, 'macs2/anno/{pair}.peakAnno.{ext}'), pair=PAIRS, ext=['pdf', 'txt']),
        # bw data quality
        expand(opj(OUTDIR, "macs2/matrix/{control}.{pos}.{types}.pdf"), control=CONTROL.keys(), pos=['tss', 'genebody'], types=['heatmap', 'profile']),
        # correlation heatmap
        opj(OUTDIR, "macs2/bigwig/heatmap_spearman_corr_readCounts.pdf"),
        # get stats
        # "report/stats.csv",