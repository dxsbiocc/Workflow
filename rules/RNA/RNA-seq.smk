import os
import json

# configfile
configfile: os.path.join(PATH, "config/RNA-seq.yaml")
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
# rsem parameters
RSEM_INDEX = config['parameters']['rsem']['index']
# salmon parameters
SALMON_INDEX = config['parameters']['salmon']['index']
# kilisto parameters
KALLISTO_INDEX = config['parameters']['kallisto']['index']
############################################################
#                          Include                         #
############################################################
include: os.path.join(PATH, "rules/common/trimmed.smk")
include: os.path.join(PATH, "rules/common/mapping.smk")
include: os.path.join(PATH, "rules/RNA/fusion.smk")
# inter-groups differential analysis
if GROUPS:
    include: os.path.join(PATH, "rules/RNA/differential_expression.smk")
############################################################
#                           Runing                         #
############################################################
rule all:
    input:
        expand(opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam"), sample=SAMPLES),