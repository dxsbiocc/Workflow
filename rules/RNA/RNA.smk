import os
import json

# configfile
configfile: os.path.join(PATH, "config/RNA-seq.yaml")
############################################################
#                      Global Variable                     #
############################################################
# sample groups
# sample groups
if 'type' in DATA.columns:
    GROUPS = DATA.groupby('type').apply(lambda x: x.index.to_list()).to_dict()
    OUTPUT.append(opj(OUTDIR, "expression/DGEs.csv"))
    OUTPUT.append(opj(OUTDIR, "expression/volcano.pdf"))
    OUTPUT.append(opj(OUTDIR, "enrichment/enrichment.rda"))
    OUTPUT.append(opj(OUTDIR, "expression/PCA.pdf"))
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
# mapping in rsem 
if QUANTIFY_TOOL != 'rsem' and config["parameters"]["rsem"]["mode"] != 'fastq':
    OUTPUT.append(expand(opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam"), sample=SAMPLES))
# fusion
FUSION =  config['control']['fusion']
if isinstance(FUSION, str):
    if FUSION.lower() == "arriba":
        OUTPUT.append(expand(opj(OUTDIR, "fusion/arriba/{sample}.tsv"), sample=SAMPLES))
    elif FUSION.lower() == "star-fusion":
        OUTPUT.append(expand(opj(OUTDIR, "fusion/star-fusion/{sample}"), sample=SAMPLES))
## star fusion
GENOME_LIB = config['parameters']['star-fusion']['genome_lib_dir']  # https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/
## arriba
BLACKLIST = config['parameters']['arriba'].get('blacklist', None)
KNOWN_FUSIONS = config['parameters']['arriba'].get('known_fusions', None)
SV_FILE = config['parameters']['arriba'].get('sv', None)
# enrichment
DBMAP = {
    'kegg': opj(PATH, 'data/pathway/kegg.pathway.csv'),
    'hsa': opj(PATH, 'data/pathway/c2.cp.v2023.1.Hs.symbols.gmt'),
    'mmu': opj(PATH, 'data/pathway/m2.cp.v2023.2.Mm.symbols.gmt')
}
ORGDB = config['parameters']['enrichment'].get('orgdb')
pathway = config['parameters']['enrichment'].get('pathway')
if os.path.isfile(pathway):
    PATHWAYDB = pathway
else:
    PATHWAYDB = DBMAP.get(pathway.lower(), '')
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
include: os.path.join(PATH, "rules/RNA/enrichment.smk")
############################################################
#                           Runing                         #
############################################################
localrules:
    mapped_stats, mapped_plotBamStats, mapped_idxstats, mapped_flagstat, \
    dedup_stats, dedup_plotBamStats, dedup_idxstats, dedup_flagstat

rule use_all:
    input:
        expand(opj(OUTDIR, "quantity/{sample}/{sample}.quant"), sample=SAMPLES),
        opj(OUTDIR, "expression/expression.tsv"),
        
    message:
        print('RNA-seq pipeline is runing!')