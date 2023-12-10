import os
import json
import getpass
import pandas as pd
from snakemake.utils import min_version
# snakemake built-in report function requires min version 7.12
min_version("7.12.0")

# configfile
configfile: os.path.join(PATH, "config/Assembly.yaml")
############################################################
#                      Global Variable                     #
############################################################
DATA = pd.read_csv(config['data']['sample_file'], index_col=0)
DATA_DICT = {name: group.index.to_list() for name, group in DATA.groupby('type')}
SAMPLES = DATA.query('type != "HIFI"').index.to_list()

ASSESS = ['HIFI']
if DATA_DICT['DNA']:
    ASSESS.append('DNA')

# assemble mode: solo/hic/trio
ASSEMBLE_MODE = config['control']['assembly_mode'].lower()
if DATA_DICT['HIC']:
    ASSEMBLE_MODE = "hic"
TRIO = []
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
    suffix = "|".join(['p_ctg', 'hap1.p_ctg', 'hap2.p_ctg'])
    hap = "|".join(['hap1.p_ctg', 'hap2.p_ctg'])

############################################################
#                          Include                         #
############################################################
include: os.path.join(PATH, "rules/Assembly/preprocessing.smk")
include: os.path.join(PATH, "rules/Assembly/assessing.smk")
include: os.path.join(PATH, "rules/Assembly/assemble.smk")
include: os.path.join(PATH, "rules/Assembly/quality_control.smk")
############################################################
#                           Runing                         #
############################################################
rule use_all:
    input:
        expand(opj(OUTDIR, "trimmed/{hifi}/"), hifi=DATA_DICT['HIFI']),
        expand(opj(OUTDIR, "trimmed/{sample}/{sample}.clean.{run}.fq.gz"), sample=SAMPLES, run=RUN),
        expand(opj(OUTDIR, "genomescope/{access}/{version}"), access=ASSESS, version=['v1', 'v2']),
        expand(opj(OUTDIR, "assembly/genome.{suffix}.gfa"), suffix=["p_ctg", "hap1.p_ctg", "hap2.p_ctg"]),
        expand(opj(OUTDIR, "genome/genome.{suffix}.stat"), suffix=["p_ctg", "hap1.p_ctg", "hap2.p_ctg"]),
        expand(opj(OUTDIR, "busco/{suffix}/"), suffix=['p_ctg', 'hap1.p_ctg', 'hap2.p_ctg']),
        expand(opj(OUTDIR, "genome/genome.{suffix}.fasta.fai"), suffix=['p_ctg', 'hap1.p_ctg', 'hap2.p_ctg']),
        expand(opj(OUTDIR, "qc/stats/{suffix}.tsv"), suffix=['p_ctg', 'hap1.p_ctg', 'hap2.p_ctg'])
        opj(OUTDIR, "merqury/diploid"),
        expand(opj(OUTDIR, "quast/{suffix}"), suffix=['p_ctg', 'hap1.p_ctg', 'hap2.p_ctg'])
        expand(opj(OUTDIR, "scaffolding/juicer/{hap}/aligned/merged_nodups.txt"), hap=['hap1.p_ctg', 'hap2.p_ctg'])
        expand(opj(OUTDIR, "scaffolding/3d-dna/{hap}"), hap=['hap1.p_ctg']),
