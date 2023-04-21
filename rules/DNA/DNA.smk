


# include utils config
include: os.path.join(PATH, "rules/common/utils.smk")

# known sites
KNOWN_SNP = ["dbsnp_146.hg38.vcf.gz"]
KNOWN_INDEL = ["Homo_sapiens_assembly38.known_indels.vcf.gz", "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"]
KNOWN = KNOWN_SNP + KNOWN_INDEL

# 1. trimming
get_trimmed(config['control']['trimming'])
# 2. mapping
get_mapping(config['control']['mapping'])
# 3. rmdup
get_dedup(config['control']['dedup'])
# index
rule dedup_index:
    input: 
        "dedup/{sample}.rmdup.bam"
    output:
        "dedup/{sample}.rmdup.bam.bai"
    log:
        "logs/dedup/dmdup_index_{sample}.log"
    wrapper:
        get_wrapper('samtools', 'index')

# 4. BQSR
include: "gatk4.smk"
# 5. call mutation
# 6.
rule combinegvcfs:
    input:
        gvcfs = expand(
            "called/{sample}.g.vcf.gz", sample=samples.index
        ),
        ref = REFERENCE,
    output:
        gvcf = "called/all.g.vcf.gz",
    log:
        "logs/called/combinegvcfs.log",
    wrapper:
        "0.74.0/bio/gatk/combinegvcfs"

rule genotypegvcfs:
    input:
        gvcf = rules.combinegvcfs.output.gvcf,
        ref = REFERENCE,
    output:
        vcf = temp("genotyped/all.vcf.gz"),
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"],
    log:
        "logs/gatk/genotypegvcfs.{contig}.log",
    wrapper:
        "0.74.0/bio/gatk/genotypegvcfs"

rule mergevcfs:
    input:
        vcfs=lambda w: expand(
            "genotyped/all.{contig}.vcf.gz", contig=get_contigs()
        ),
    output:
        vcf="results/genotyped/all.vcf.gz",
    log:
        "logs/picard/merge-genotyped.log",
    wrapper:
        "0.74.0/bio/picard/mergevcfs"
