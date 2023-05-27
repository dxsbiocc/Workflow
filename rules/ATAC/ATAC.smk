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
if 'type' in samples.columns:
    MACS2_MAP = samples['type'].to_dict()

include: os.path.join(PATH, "rules/common/trimmed.smk")
include: os.path.join(PATH, "rules/ATAC/common.smk")
include: os.path.join(PATH, "rules/ATAC/mapping.smk")
include: os.path.join(PATH, "rules/ATAC/filtering.smk")
include: os.path.join(PATH, "rules/ATAC/callpeak.smk")
include: os.path.join(PATH, "rules/ATAC/motif.smk")
include: os.path.join(PATH, "rules/ATAC/quality.smk")
include: os.path.join(PATH, "rules/ATAC/report.smk")

# ------------------------- run rules ---------------------- #
# report: "report/workflow.rst"
# local
localrules:
    stats, idxstats, flagstat

rule use_all:
    input:
        # trimmed adapter
        expand("trimmed/{sample}/{sample}.clean.{unit}.fq.gz", sample=SAMPLES, unit=['R1', 'R2']),
        # spikein
        expand("bowtie2/{sample}/{sample}.spikein.bam", sample=SAMPLES),
        # stats
        expand("bowtie2/stats/{sample}.{stats}", sample=SAMPLES, stats=['stats', 'idxstats', 'flagstats']),
        expand("bowtie2/plot/{sample}", sample=SAMPLES),
        # mapping and remove duplicates
        expand("dedup/{sample}/{sample}.cons.bam.bai", sample=SAMPLES),
        # shift 9bp
        expand("dedup/{sample}/{sample}.filtered.bam", sample=SAMPLES),
        # callpeak
        expand("macs2/narrow/{pair}_peaks.narrowPeak", pair=PAIRS),
        expand("macs2/broad/{pair}_peaks.broadPeak", pair=PAIRS),
        # convert to bigwig
        expand("macs2/bigwig/{sample}.norm.bw", sample=SAMPLES),
        # extract peaks
        expand("macs2/narrow/{pair}_homer_peaks.txt", pair=PAIRS),
        # motify analysis
        expand("motifs/{pair}", pair=PAIRS),
        # peak annotate
        expand('macs2/anno/{pair}.peakAnno.{ext}', pair=PAIRS, ext=['pdf', 'txt']),
        # bw data quality
        expand("macs2/matrix/{control}.{pos}.{types}.pdf", control=CONTROL.keys(), pos=['tss', 'genebody'], types=['heatmap', 'profile']),
        # correlation heatmap
        "macs2/bigwig/heatmap_spearman_corr_readCounts.pdf",
        # get stats
        "report/stats.csv",