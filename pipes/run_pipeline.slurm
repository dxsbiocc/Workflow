#!/bin/bash
# example run snakemake workflow
# use of cluster: SLURM, sample file, rules in modules: cf Snakefile_modules.smk
# using a config file in yaml
# using ressources CPU and memory
# Using cluster with SLURM cluster parametrer in a cluster.json
# all binaries are in singularity images
# Samples (samples names and fastq files) in csv file

#SBATCH --job-name=pipeline # job name (-J)
#SBATCH --time="72:00:00" #max run time "hh:mm:ss" or "dd-hh:mm:ss" (-t)
#SBATCH --cpus-per-task=1 # max nb of cores (-c)
#SBATCH --ntasks=1 #nb of tasks
#SBATCH --mem=64G # max memory (-m)
#SBATCH --output=pipeline.%j.out # stdout (-o)

########################## On genotoul ###############################
# uncomment the modules
## snakemake 5.3
#module load system/Python-3.6.3
#module load system/singularity-3.5.3
######################################################################

########################## On ifb-core ###############################
## uncomment the next line
#module load snakemake
######################################################################


# Using external rules in the main Snakefile (minimal Snakefile)
RULES=Snakefile
# Or with all rules included in one snakemake file:
#RULES=Snakefile.smk

# config file to be edited
CONFIG=../config/config.yaml

# slurm directive by rule (can be edited if needed)
CLUSTER_CONFIG=cluster.yaml
# sbatch directive to pass to snakemake
CLUSTER='sbatch --mem={cluster.mem} -t {cluster.time} -c {cluster.cores} -J {cluster.jobname} --nodelist={cluster.nodelist} -o logs/{cluster.out} -e logs/{cluster.error}'
# Maximum number of jobs to be submitted at a time (see cluster limitation)
MAX_JOBS=500

# Full clean up:
#rm -fr .snakemake logs *.out *.html *.log *.pdf
# Clean up only the .snakemake
rm -fr .snakemake

# log directory is mandatory (see $CLUSTER) else slurm jobs failed but not the master job
mkdir -p logs

# Generate the dag files
# With samples
snakemake --configfile $CONFIG -s $RULES --dag | dot -Tpdf > dag.pdf
# Only rules
# snakemake --configfile $CONFIG -s $RULES --rulegraph | dot -Tpdf > dag_rules.pdf

# Dry run (simulation)
# snakemake --configfile $CONFIG -s $RULES -np -j $MAX_JOBS --cluster-config $CLUSTER_CONFIG --cluster "$CLUSTER" >snakemake_dryrun.out

# Full run (if everething is ok: uncomment it)
snakemake --configfile $CONFIG -s $RULES -p -j $MAX_JOBS --cluster-config $CLUSTER_CONFIG --cluster "$CLUSTER"

# If latency problem add to the run:
# --latency-wait 60

exit 0

