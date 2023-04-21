from collections import defaultdict

if SAMPLE_MAP:
    CONTROL = defaultdict(list)
    for k, v in pairs.items():
        CONTROL[v].append(k)
else:
    CONTROL = {k : '', for k in PAIRS}


include: os.path.join(PATH, "rules/common/utils.smk")
# get rules
get_trimmed(config['control']['trimmed_tool'])
# get_mapping(config['control']['mapping_tool'])

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

rule all:
    input:
        # trimmed adapter
        expand("trimmed/{sample}.clean.{unit}.fq.gz", sample=SAMPLES, unit=['R1', 'R2']),
        # spikein
        expand("bowtie2/{sample}.spikein.bam", sample=SAMPLES),
        # stats
        expand("bowtie2/stats/{sample}.{stats}", sample=SAMPLES, stats=['stats', 'idxstats', 'flagstats']),
        expand("bowtie2/plot/{sample}", sample=SAMPLES),
        # mapping and remove duplicates
        expand("dedup/{sample}.cons.bam.bai", sample=SAMPLES),
        # shift 9bp
        expand("dedup/{sample}.filtered.bam", sample=SAMPLES),
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
        # get stats
        "report/stats.csv",