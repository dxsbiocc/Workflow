# include utils config
include: os.path.join(PATH, "rules/common/utils.smk")

DIGESTION = 'single'  # or 'double'

# Snakefile

# Configuration
configfile: "config.yaml"

# Define input files
reads = config["reads"]
ref_genome = config["ref_genome"]
fragments = config["fragments"]

# 1. trimming
get_trimmed(config['control']['trimming'])
# 2. mapping
get_mapping(config['control']['mapping'])
# 3. rmdup
get_dedup(config['control']['dedup'])
# 4. HiCExplorer
include: os.path.join(PATH, "rules/HiC/hicexplorer.smk")
# 5. multiQC
include: os.path.join(PATH, "rules/common/multiqc.smk")

onstart:
    if config["verbose"]:
        print("--- Workflow parameters --------------------------------------------------------")
        print("samples:", SAMPLES)
        print("genome index:", INDEX)
        print("-" * 80, "\n")

        print("--- Environment ----------------------------------------------------------------")
        print("$TMPDIR: ",os.getenv('TMPDIR', ""))
        print("$HOSTNAME: ",os.getenv('HOSTNAME', ""))
        print("-" * 80, "\n")

# Rules
rule all:
    input:
        "results/filtered_reads.fastq",

onsuccess:
    if config["verbose"]:
        print("\n--- Hi-C workflow finished successfully! --------------------------------\n")

onerror:
    print("\n !!! ERROR in HI-C workflow! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")