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
else:
    GROUPS = None
# differential analysis metohd
DIFF_TOOL = config['control']['differential_analysis'].lower()
# star parameters
READ_LENGTH = config['parameters']['star']['read_length']
GENOME_LIB = config['parameters']['star']['genome_lib_dir']
# quantify tool
QUANTIFY_TOOL = config['control']['quantify'].lower()
# quantify index
QUANTIFY_INDEX = config['parameters'][QUANTIFY_TOOL].get("index", None)
############################################################
#                          Include                         #
############################################################
include: os.path.join(PATH, "rules/common/trimmed.smk")
include: os.path.join(PATH, "rules/common/mapping.smk")
include: os.path.join(PATH, "rules/common/stats.smk")
include: os.path.join(PATH, "rules/common/dedup.smk")
include: os.path.join(PATH, "rules/RNA/quantify.smk")
include: os.path.join(PATH, "rules/RNA/fusion.smk")
# inter-groups differential analysis
if GROUPS:
    include: os.path.join(PATH, "rules/RNA/expression.smk")
############################################################
#                           Runing                         #
############################################################
rule use_all:
    input:
        expand(opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam"), sample=SAMPLES),
        expand(opj(OUTDIR, "quantity/{sample}"), sample=SAMPLES),
    message:
        print('RNA-seq pipeline is runing!')