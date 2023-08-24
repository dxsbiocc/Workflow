import os
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version


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

def get_environment(*args, local=True):
    """Get environment.yaml"""
    if local:
        abspath = os.path.join(PATH, 'wrappers', *args, 'environment.yaml')
        return abspath
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
        
# ------------------------------ DNA ------------------------------ #
def get_filter(mode):
    if mode == 'snp':
        filters = {"my_filter": SNP_FILTER}
    elif mode == 'indel':
        filters = {"my_filter": INDEL_FILTER}
    else:
        raise ValueError('Unsupport filter!')
    return filters
    
def get_contig():
    if 'REGION_BED' in config:
        pass
    else:
        with open(f'{REFERENCE}.fai') as fp:
            contig = pd.read_csv(fp, sep='\t', header=None, usecols=[0], dtype=str).squeeze().to_list()   
        return contig