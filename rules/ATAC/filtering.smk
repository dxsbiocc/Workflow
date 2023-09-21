if SHIFT:

    rule shift:
        input:
            rules.mark_duplicates.output.bam,
            rules.cons_index.output,
        output:
            "dedup/{sample}/{sample}.shift.bam",
        log:
            "logs/dedup/dedup_shift_{sample}.log",
        threads: 4
        params:
            # optional parameters
            extra="{shift} --samFlagInclude {flag} --blackListFileName {bl} --minMappingQuality {mapq}".format(
                shift=config["parameters"]["alignmentSieve"]["shift"],
                flag=config["parameters"]["alignmentSieve"]["flag"],
                bl=BLACKLIST,
                mapq=config["parameters"]["alignmentSieve"]["mapq"],
            ),
        wrapper:
            get_wrapper("deeptools", "alignmentsieve")

else:

    rule noshift:
        input:
            rules.mark_duplicates.output.bam,
            rules.cons_index.output,
        output:
            "dedup/{sample}/{sample}.shift.bam",
        log:
            "logs/dedup/dedup_no_shift_{sample}.log",
        threads: 4
        params:
            blacklist=BLACKLIST,
        conda:
            lambda wildcards: get_environment("bedtools", "intersect")
        shell:
            "bedtools intersect -v -abam {input[0]} -b {params.blacklist} > {output}"


rule shift_sort:
    input:
        "dedup/{sample}/{sample}.shift.bam",
    output:
        "dedup/{sample}/{sample}.shift.sort.bam",
    log:
        "logs/dedup/dedup_shift_sort_{sample}.log",
    params:
        extra="-m 4G",
    threads: 8
    wrapper:
        get_wrapper("samtools", "sort")


rule filterChrM:
    input:
        rules.shift_sort.output,
    output:
        "dedup/{sample}/{sample}.filtered.bam",
    conda:
        lambda wildcards: get_environment("samtools", "view")
    shell:
        "(samtools view -h {input} | grep -v chrM | samtools view -bS -F 0x4 -q 30 | samtools sort -@ 4 -o {output}) "
        "&& (samtools index {output})"
