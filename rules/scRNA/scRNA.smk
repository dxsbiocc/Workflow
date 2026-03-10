import os
import json


############################################################
#                      Global Variable                     #
############################################################
# aligner
assert MAPPING in ['star', 'cellranger', 'simpleaf', 'kb-python'], \
                    'aligner must be one of starsolo, cellranger, simpleaf, kb-python'
if MAPPING == 'star':
    ALIGNER = 'starsolo'
else:
    ALIGNER = MAPPING
# platform
PLATFORM = config[ALIGNER]['platform']
# star parameters
READ_LENGTH = config['parameters']['starsolo']['read_length']
# cellranger parameters
CELLRANGER = config['parameters']['cellranger']['executable']
# simpleaf parameters

# kb-python parameters

# add rules
OUTPUT.append(expand(opj(OUTDIR, ALIGNER, "{sample}"), sample=SAMPLES))
############################################################
#                          Include                         #
############################################################
if TRIMMING:
    include: os.path.join(PATH, "rules/common/trimmed.smk")

include: os.path.join(PATH, f"rules/scRNA/{ALIGNER}.smk")
############################################################
#                           Runing                         #
############################################################

rule use_all:
    input:
        expand(opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam"), sample=SAMPLES),
    message:
        print('RNA-seq pipeline is runing!')