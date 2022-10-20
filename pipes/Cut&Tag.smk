# Snakemake workflow example.
#
# Snakemake feature: fastq from csv file, config file in yaml, modules (rules), SLURM (json), report


import pandas as pd
from snakemake.utils import min_version
# snakemake built-in report function requires min version 5.1
min_version("5.1.0")

configfile: "config.yaml"

#read the sample file using pandas lib (sample names+ fastq names) and crezate index using the sample name
samples_info = pd.read_csv(config["samples_info"], index_col=0, dtype=str, comment='#')

#List unique values in the samples['Group'] column
groups = samples.Group.unique() #list of uniq group name


##### load rules #####
include: "rules/utils.smk"
# work dir
workdir: config['workdir']


rule all:
    input:

        expand("{outdir}/multiqc/multiqc_report.html", outdir=config["outdir"]),
        expand("{outdir}/mapped/{sample}_sorted.bam", outdir=config["outdir"], sample=samples['SampleName']),
        "{outdir}/bams.list".format(outdir=config["outdir"])
        #expand("{outdir}/mapped/{group}_sorted.bam", outdir=config["outdir"], group=groups)
