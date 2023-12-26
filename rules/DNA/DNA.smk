import os

# configfile
# configfile: os.path.join(PATH, "config/DNA-seq.yaml")
############################################################
#                      Global Variable                     #
############################################################
# control parameters
MODE = ['snp', 'indel']
FILTER = config['control']['filter']
# --------------------- Reference data ------------------- #
DICT = config['data']['dict']
# ------------------------- Region ----------------------- #
REGION_BED = config['data']['REGION_BED'] if 'REGION_BED' in config['data'] else ''

CONTIG = get_contig()
# ---------------------- Known sites --------------------- #
KNOWEN_SITE = [config['data']["known_site"][k] for k in config['data']["known_site"]]
KNOWEN_SITE_DBSNP = config['data']["known_site"]['dbsnp']
# ---------------------- VQSR filters -------------------- #
MILLS = config['parameters']['gatk']['vqsr']["mills"]
MILLS_IDX = config['parameters']['gatk']['vqsr']["mills_idx"]
OMNI = config['parameters']['gatk']['vqsr']["omni"]
OMIN_IDX = config['parameters']['gatk']['vqsr']["omni_idx"]
G1K = config['parameters']['gatk']['vqsr']["g1k"]
G1K_IDX = config['parameters']['gatk']['vqsr']["g1k_idx"]
DBSNP = config['parameters']['gatk']['vqsr']["dbsnp"]
DBSNP_IDX = config['parameters']['gatk']['vqsr']["dbsnp_idx"]
# ---------------------- Hard filters -------------------- #
SNP_FILTER = config['parameters']['gatk']['hard_filter']['snp_filter']
INDEL_FILTER = config['parameters']['gatk']['hard_filter']['indel_filter']
# --------------------- VEP parameters ------------------- #
VEP_CACHE = config['parameters']['vep']['local']['cache']
VEP_PLUGINS_LOCAL = config['parameters']['vep']['local']['plugins']
VEP_SPECIES = config['parameters']['vep']['online']['species']
VEP_BUILD = config['parameters']['vep']['online']['build']
VEP_RELEASE = config['parameters']['vep']['online']['release']
VEP_PLUGINS = config['parameters']["vep"]["plugins"]
# -------------------- Delly parameters ------------------ #
DELLY_EXCLUDE = config['parameters']['delly']['exclude']
# ------------------ Wildcard constraints ---------------- #
wildcard_constraints:
    contig = "|".join(CONTIG),
    mode = "snp|indel"

############################################################
#                         0. Include                       #
############################################################
include: os.path.join(PATH, "rules/DNA/check.smk")
include: os.path.join(PATH, "rules/common/trimmed.smk")
include: os.path.join(PATH, "rules/common/mapping.smk")
include: os.path.join(PATH, "rules/common/dedup.smk")
include: os.path.join(PATH, "rules/common/stats.smk")
include: os.path.join(PATH, "rules/DNA/gatk4.smk")
include: os.path.join(PATH, "rules/DNA/filtering.smk")
include: os.path.join(PATH, "rules/DNA/anno.smk")
include: os.path.join(PATH, "rules/DNA/sv.smk")
include: os.path.join(PATH, "rules/DNA/cnv.smk")
include: os.path.join(PATH, "rules/DNA/report.smk")
############################################################
#                           Runing                         #
############################################################
rule use_all:
    input:
        # data process
        expand(opj(OUTDIR, "trimmed/{sample}/{sample}.clean.{run}.fq.gz"), sample=SAMPLES, run=RUN),
        expand(opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam"), sample=SAMPLES),
        expand(opj(OUTDIR, "dedup/{sample}/{sample}.rmdup.bam"), sample=SAMPLES),
        expand(opj(OUTDIR, "report/{bam}/stats/{sample}.{stats}"), bam=["mapped", "dedup"], sample=SAMPLES, stats=['stats', 'idxstats', 'flagstats']),
        expand(opj(OUTDIR, "report/{bam}/plot/{sample}"), bam=["mapped", "dedup"], sample=SAMPLES),
        # variant called
        expand(opj(OUTDIR, "variant/filtered/all.{mode}.filtered.vcf.gz"), mode=MODE),
        # annotate
        opj(OUTDIR, "annotated/all.vcf.gz"),
        opj(OUTDIR, "sv/all.vcf.gz")