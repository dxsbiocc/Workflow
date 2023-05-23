import os
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version


# -------------------- Global Parameters -------------------- #
# database path
if config['data']['db'] == 'local':
    DATABASE = os.path.join(PATH, 'data')
else:
    DATABASE = config['data']['db']
# paired
if config['control']['paired']:
    RUN = ['R1', 'R2']
else:
    RUN = ['R1']

# reference
REFERENCE = config['data']['ref']
INDEX = config['data']['index']
GTF = config['data']['gtf']
MAPPING = config['control']['mapping']
TRIMMING = config['control']['trimming']
DEDUP = config['control']['dedup']
# ------------------------- modules ------------------------ #
# define modules
module trimmed_workflow:
    snakefile:
        os.path.join(PATH, "rules/common/trimmed.smk")

module mapping_workflow:
    snakefile:
        os.path.join(PATH, "rules/common/mapping.smk")

module dedup_workflow:
    snakefile:
        os.path.join(PATH, "rules/common/dedup.smk")
# import rules
def get_trimmed(rule_name):
    if rule_name == "fastp":
        use rule fastp from trimmed_workflow
    elif rule_name == "trimmomatic":
        use rule trimmomatic from trimmed_workflow
    elif rule_name == "cutadapt":
        use rule cutadapt from trimmed_workflow
    elif rule_name == "trim_galore":
        use rule trim_galore from trimmed_workflow
    else:
        raise ValueError(f'the rule: {rule_name} not support!')

def get_mapping(rule_name):
    if rule_name == "bwa":
        use rule bwa from mapping_workflow
    elif rule_name == "bowtie2":
        use rule bowtie2 from mapping_workflow
    elif rule_name == "hisat2":
        use rule hisat2 from mapping_workflow
    elif rule_name == "star":
        use rule star from mapping_workflow
    elif rule_name == "minimap2":
        use rule minimap2 from mapping_workflow
    else:
        raise ValueError(f'the rule: {rule_name} not support!')

def get_dedup(rule_name):
    if rule_name == "markduplicates":
        use rule markduplicates from dedup_workflow
    elif rule_name == "sambamba":
        use rule sambamba from dedup_workflow
    else:
        raise ValueError(f'the rule: {rule_name} not support!')
# -------------------- Global Functions -------------------- #
def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = samples.loc[wildcards.sample, ["fastq1", "fastq2"]].dropna()
    if config['control']['paired']:
        return [fastqs.fastq1, fastqs.fastq2]
    return [fastqs.fastq1]

def get_wrapper(*args, local=True):
    """Get wrappers path"""
    if local:
        abspath = os.path.join(PATH, 'wrappers', *args)
        return "file://{}".format(abspath)
    else:
        raise ValueError("Please use local version!")

def get_script(script):
    return os.path.join(PATH, 'scripts', script)

def get_adapter(method='fastp'):
    """get adapter path"""
    if config['control']['adapters'] == 'Truseq':
        adapt_file = os.path.join(DATABASE, "adapters", "TruSeqAdapters.fa")
    elif config['control']['adapters'] == 'Nextera':
        adapt_file = os.path.join(DATABASE, "adapters", "NexteraPE-PE.fa")
    else:
        raise ValueError('%s not support!' % config['control']['adapters'])

    if method == 'fastp':
        pat = f"--adapter_fasta {adapt_file}"
    elif method == 'trimmomatic':
        pat = f'ILLUMINACLIP:{adapt_file}:2:15:4:4:true'
    elif method == 'cutadapt':
        pat = f'-a "file:{adapt_file}"'
    else:
        raise ValueError(f'{method} not support!')
    return pat