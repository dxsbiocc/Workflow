import os
from snakemake.logging import logger

############################################################
#                      Global Variable                     #
############################################################
if 'group' in DATA.columns:
    GROUPS = DATA.groupby('group').apply(lambda x: x.index.to_list()).to_dict()
    GROUPS_LINE = '\t'.join(DATA['group'].values)
else:
    GROUPS = None
    GROUPS_LINE = ""

# ----------------- reads ------------------
REMOVE_HOST = config['control'].get('remove_host')
READS = config['parameters']['humann'].get('reads', None)
TAXONOMY = config['parameters']['bracken'].get('taxonomy', ['S'])
# ---------------- assembly ----------------
assembly_method = config['control'].get('assembly_method', 'megahit')
if assembly_method == "metaspades":
    if len(SAMPLES) > 1:
        logger.warning("MetaSPAdes is designed primarily for single-sample metagenomic assembly. For multiple samples, we recommend using megahit.")
    OUTPUT.append(expand(opj(OUTDIR, "metaspades/{sample}/final.contigs.fasta"), sample=SAMPLES))
    ASSEMBLY = ["metaspades"]
elif assembly_method == "megahit":
    OUTPUT.append(expand(opj(OUTDIR, "megahit/{sample}/final.contigs.fasta"), sample=SAMPLES))
    ASSEMBLY = ["megahit"]
elif assembly_method == "both":
    if len(SAMPLES) > 1:
        logger.warning("MetaSPAdes is designed primarily for single-sample metagenomic assembly. For multiple samples, we recommend using megahit.")
    OUTPUT.append(expand(opj(OUTDIR, "megahit/{sample}/final.contigs.fasta"), sample=SAMPLES))
    OUTPUT.append(expand(opj(OUTDIR, "metaspades/{sample}/final.contigs.fasta"), sample=SAMPLES))
    OUTPUT.append(expand(opj(OUTDIR, "quast/compare/{sample}/report.html"), sample=SAMPLES))
    ASSEMBLY = ["megahit", "metaspades"]
else:
    raise ValueError(f"Invalid assembly method: {assembly_method}")

quast_method = config['control'].get('quast', 'quast')
if quast_method == "quast":
    OUTPUT.append(expand(opj(OUTDIR, "quast/{assembly}/{sample}/report.html"), assembly=ASSEMBLY, sample=SAMPLES))
elif quast_method == "metaquast":
    OUTPUT.append(expand(opj(OUTDIR, "metaquast/{assembly}/{sample}/report.html"), assembly=ASSEMBLY, sample=SAMPLES))
else:
    raise ValueError(f"Invalid quast method: {quast_method}")

FULL_GENE = config['control'].get('full_gene', True)
# ---------------- binning ----------------
BINNING_WITHOUT_ASSEMBLY = config['control'].get('binning_without_assembly', True)
if config['parameters']['metawrap']['assembly_method'] not in ASSEMBLY:
    BINNING_WITHOUT_ASSEMBLY = False
############################################################
#                         0. Include                       #
############################################################
include: os.path.join(PATH, "rules/common/trimmed.smk")
include: os.path.join(PATH, "rules/Meta/filtering.smk")
include: os.path.join(PATH, "rules/Meta/reads.smk")
include: os.path.join(PATH, "rules/Meta/assembly.smk")
include: os.path.join(PATH, "rules/Meta/binning.smk")
############################################################
#                           Runing                         #
############################################################
localrules:
    humann_associate, humann_barplot, humann_regroup_ko, humann_join_ko,
    humann_join, humann_split_abundance, graphlan_export, graphlan_annotate,
    graphlan_graphlan

rule use_all:
    input:
        expand(opj(OUTDIR, "trimmed/{sample}/{sample}.clean.{unit}.fq.gz"), sample=SAMPLES, unit=RUN),
        expand(opj(OUTDIR, "downsample/{sample}.sampled.fq.gz"), sample=SAMPLES),
        expand(opj(OUTDIR, "humann/{sample}"), sample=SAMPLES),
        # humann
        opj(OUTDIR, "humann/metaphlan4/taxonomy.tsv"),
        opj(OUTDIR, "humann/abundance.tsv"),
        opj(OUTDIR, "humann/abundance.relative.tsv"),
        opj(OUTDIR, "humann/abundance.stratified.tsv"),
        opj(OUTDIR, "humann/abundance.unstratified.tsv"),
        opj(OUTDIR, "humann/abundance.pcl"),
        opj(OUTDIR, "humann/associate.tsv"),
        opj(OUTDIR, "humann/barplot"),
        expand(opj(OUTDIR, "humann/ko/{sample}.ko.tsv"), sample=SAMPLES),
        opj(OUTDIR, "humann/ko/ko.tsv"),
        opj(OUTDIR, "humann/ko/ko.stratified.tsv"),
        opj(OUTDIR, "humann/ko/ko.unstratified.tsv"),
        # graphlan
        opj(OUTDIR, "graphlan/graphlan.pdf"),
        # lefse
        opj(OUTDIR, "lefse/cladogram.pdf"),
        opj(OUTDIR, "lefse/result.pdf"),
        opj(OUTDIR, "lefse/features"),
        # kraken2
        expand(opj(OUTDIR, "kraken2/{sample}/{sample}.report"), sample=SAMPLES),
        expand(opj(OUTDIR, "kraken2/{sample}/{sample}.mpa"), sample=SAMPLES),
        expand(opj(OUTDIR, "kraken2/{sample}/{sample}.percentages.mpa"), sample=SAMPLES),
        opj(OUTDIR, "kraken2/merged.mpa.tsv"),
        # bracken
        expand(opj(OUTDIR, "bracken/{taxonomy}/{sample}.bracken"), taxonomy=TAXONOMY, sample=SAMPLES),
        expand(opj(OUTDIR, "bracken/{taxonomy}/merged.bracken.tsv"), taxonomy=TAXONOMY),
        expand(opj(OUTDIR, "bracken/{taxonomy}/alpha.tsv"), taxonomy=TAXONOMY),
        expand(opj(OUTDIR, "bracken/{taxonomy}/beta.txt"), taxonomy=TAXONOMY),
        expand(opj(OUTDIR, "bracken/{taxonomy}/{sample}.krona.html"), taxonomy=TAXONOMY, sample=SAMPLES),
        # NR
        expand(opj(OUTDIR, "NR/{assembly}/nucleotide.fa"), assembly=ASSEMBLY),
        expand(opj(OUTDIR, "NR/{assembly}/protein.fa"), assembly=ASSEMBLY),
        expand(opj(OUTDIR, "NR/{assembly}/protein_no_end.fa"), assembly=ASSEMBLY),
        expand(opj(OUTDIR, "NR/{assembly}/all_gene.fa"), assembly=ASSEMBLY),
        expand(opj(OUTDIR, "NR/{assembly}/gene.fa"), assembly=ASSEMBLY),
        # salmon
        expand(opj(OUTDIR, "salmon/{assembly}/quant.tpm"), assembly=ASSEMBLY),
        expand(opj(OUTDIR, "salmon/{assembly}/quant.counts"), assembly=ASSEMBLY),
        # eggnog
        expand(opj(OUTDIR, "eggnog/{assembly}/eggnog.annotations"), assembly=ASSEMBLY),
        expand(opj(OUTDIR, "rgi/{assembly}/protein.json"), assembly=ASSEMBLY),
        expand(opj(OUTDIR, "rgi/{assembly}/nucleotide.json"), assembly=ASSEMBLY),
        # metawrap
        expand(opj(OUTDIR, "metawrap/binning/{sample}"), sample=SAMPLES),
        expand(opj(OUTDIR, "metawrap/bins_refine/{sample}"), sample=SAMPLES),
        # output files
        OUTPUT
    message:
        print('Meta pipeline is runing!')