import os
import json

# configfile
# configfile: os.path.join(PATH, "config/RNA-seq.yaml")
############################################################
#                      Global Variable                     #
############################################################
# sample groups
if 'type' in DATA.columns:
    GROUPS = DATA.groupby('type').apply(lambda x: x.index.to_list()).to_dict()
    OUTPUT.append(opj(OUTDIR, "expression/DGEs.csv"))
    OUTPUT.append(opj(OUTDIR, "expression/volcano.pdf"))
else:
    GROUPS = None
# differential analysis metohd
DIFF_TOOL = config['control']['differential_analysis'].lower()
# star parameters
READ_LENGTH = config['parameters']['star']['read_length']
# quantify tool
QUANTIFY_TOOL = config['control']['quantify'].lower()
# quantify index
if QUANTIFY_TOOL in ['salmon', 'kallisto']:
    QUANTIFY_INDEX = config['data'].get('rna_index', None)
else:
    QUANTIFY_INDEX = config['parameters'][QUANTIFY_TOOL].get("index", None)
# fusion
FUSION =  config['control']['fusion'].lower()
if FUSION == "arriba":
    OUTPUT.append(expand(opj(OUTDIR, "fusion/arriba/{sample}.tsv"), sample=SAMPLES))
elif FUSION == "star-fusion":
    OUTPUT.append(expand(opj(OUTDIR, "fusion/star-fusion/{sample}"), sample=SAMPLES))
## star fusion
GENOME_LIB = config['parameters']['star-fusion']['genome_lib_dir']  # https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/
## arriba
BLACKLIST = config['parameters']['arriba'].get('blacklist', None)
KNOWN_FUSIONS = config['parameters']['arriba'].get('known_fusions', None)
SV_FILE = config['parameters']['arriba'].get('sv', None)
############################################################
#                          Include                         #
############################################################
include: os.path.join(PATH, "rules/common/trimmed.smk")
include: os.path.join(PATH, "rules/common/mapping.smk")
include: os.path.join(PATH, "rules/common/stats.smk")
include: os.path.join(PATH, "rules/common/dedup.smk")
include: os.path.join(PATH, "rules/RNA/quantify.smk")
include: os.path.join(PATH, "rules/RNA/fusion.smk")
include: os.path.join(PATH, "rules/RNA/expression.smk")
############################################################
#                           Runing                         #
############################################################
localrules:
    mapped_stats, mapped_plotBamStats, mapped_idxstats, mapped_flagstat, \
    dedup_stats, dedup_plotBamStats, dedup_idxstats, dedup_flagstat

rule use_all:
    input:
        expand(opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam"), sample=SAMPLES),
        expand(opj(OUTDIR, "quantity/{sample}/{sample}.quant"), sample=SAMPLES),
        opj(OUTDIR, "expression/expression.tsv"),
        opj(OUTDIR, "expression/PCA.pdf"),
        OUTPUT
    message:
        print('RNA-seq pipeline is runing!')