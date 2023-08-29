import os
import json
import getpass
import pandas as pd
from snakemake.utils import min_version
# snakemake built-in report function requires min version 5.1
min_version("7.12.0")


PATH = '/data/dengxsh/WorkSpace/Projects/Workflow/'

#read the sample file using pandas lib (sample names+ fastq names) and crezate index using the sample name
DATA = pd.read_csv('/data/dengxsh/WorkSpace/assembly/sample_list.csv', index_col=0)

DATA_DICT = {name: group.index.to_list() for name, group in DATA.groupby('type')}
SAMPLES = DATA.query('type != "HIFI"').index.to_list()

ASSESS = ['HIFI']
if DATA_DICT['DNA']:
    ASSESS.append('DNA')

# 
ASSEMBLE_MODE = config['control']['assembly_mode'].lower()
if DATA_DICT['HIC']:
    ASSEMBLE_MODE = "hic"
if ASSEMBLE_MODE == 'trio':
    TRIO = DATA.query('type != "TRIO"').index.to_list()
# ------------------ Wildcard constraints ------------------ #
RUN = ['R1', 'R2']
wildcard_constraints:
    sample = "|".join(SAMPLES),
    hifi = "|".join(DATA_DICT['HIFI']),
    dna = "|".join(DATA_DICT['DNA']),
    rna = "|".join(DATA_DICT['RNA']),
    hic = "|".join(DATA_DICT['HIC']),
    access = "|".join(ASSESS),
    trio = "|".join(TRIO)

############################################################
#                         0. Include                       #
############################################################
include: os.path.join(PATH, "rules/assembly/preprocessing.smk")
include: os.path.join(PATH, "rules/assembly/assessing.smk")
include: os.path.join(PATH, "rules/assembly/assemble.smk")
############################################################
#                           Runing                         #
############################################################
rule all:
    input:
        expand("trimmed/{hifi}/", hifi=DATA_DICT['HIFI']),
        expand("trimmed/{sample}/{sample}.clean.{run}.fq.gz", sample=SAMPLES, run=RUN),
        expand("genomescope/{version}/{dna}/", version=['v1', 'v2'], dna=DATA_DICT['DNA']),
