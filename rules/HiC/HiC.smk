import os

# configfile
configfile: os.path.join(PATH, "config/HIC.yaml")
############################################################
#                      Global Variable                     #
############################################################
DIGESTION = 'single'  # or 'double'
reads = config["reads"]
ref_genome = config["ref_genome"]
fragments = config["fragments"]
############################################################
#                           Include                        #
############################################################
# 1. trimming
include: os.path.join(PATH, "rules/common/trimmed.smk")
# 2. mapping
include: os.path.join(PATH, "rules/common/mapping.smk")
# 3. rmdup
include: os.path.join(PATH, "rules/common/dedup.smk")
# 4. HiCExplorer
include: os.path.join(PATH, "rules/HiC/hicexplorer.smk")
# 5. multiQC
include: os.path.join(PATH, "rules/common/multiqc.smk")
############################################################
#                           Runing                         #
############################################################
onstart:
    if config["verbose"]:
        print("Workflow Parameters".center(80, '-'))
        print("samples:", DATA.index.to_list())
        print("genome index:", INDEX)
        print("-" * 80, "\n")

        print("Environment".center(80, '-'))
        print("$TMPDIR: ",os.getenv('TMPDIR', ""))
        print("$HOSTNAME: ",os.getenv('HOSTNAME', ""))
        print("-" * 80, "\n")

# Rules
rule use_all:
    input:
        # data process
        expand("trimmed/{sample}/{sample}.clean.{run}.fq.gz", sample=SAMPLES, run=RUN),
        expand("mapped/{sample}/{sample}.sorted.bam", sample=SAMPLES),
        expand("dedup/{sample}/{sample}.rmdup.bam", sample=SAMPLES),
        expand("report/stats/{sample}.{stats}", sample=SAMPLES, stats=['stats', 'idxstats', 'flagstats']),
        expand("report/plot/{sample}", sample=SAMPLES),

onsuccess:
    if config["verbose"]:
        print("Hi-C workflow finished successfully!".center(80, '-'))

onerror:
    print("ERROR in HI-C workflow".center(80, '!'))