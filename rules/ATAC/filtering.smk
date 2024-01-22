if SHIFT:
    rule shift:
        input:
            rules.markduplicates.output.bam,
            rules.rmdup_index.output,
        output:
            opj(OUTDIR, "dedup/{sample}/{sample}.shift.bam"),
        log:
            opj(OUTDIR, "logs/dedup/dedup_shift_{sample}.log"),
        threads: 4
        params:
            # optional parameters
            extra = "{shift} --samFlagInclude {flag} --blackListFileName {bl} --minMappingQuality {mapq}".format(
                shift=config["parameters"]['deeptools']["alignmentSieve"]["shift"],
                flag=config["parameters"]['deeptools']["alignmentSieve"]["flag"],
                bl=BLACKLIST,
                mapq=config["parameters"]['deeptools']["alignmentSieve"]["mapq"],
            ),
        wrapper:
            get_wrapper("deeptools", "alignmentsieve")
else:
    rule noshift:
        input:
            rules.markduplicates.output.bam,
            rules.rmdup_index.output,
        output:
            opj(OUTDIR, "dedup/{sample}/{sample}.shift.bam"),
        log:
            opj(OUTDIR, "logs/dedup/dedup_no_shift_{sample}.log"),
        threads: 4
        params:
            blacklist = BLACKLIST,
        conda:
            lambda wildcards: get_environment("bedtools", "intersect")
        shell:
            "bedtools intersect -v -abam {input[0]} -b {params.blacklist} > {output}"

rule shift_sort:
    input:
        opj(OUTDIR, "dedup/{sample}/{sample}.shift.bam"),
    output:
        opj(OUTDIR, "dedup/{sample}/{sample}.shift.sort.bam"),
    log:
        opj(OUTDIR, "logs/dedup/dedup_shift_sort_{sample}.log"),
    params:
        extra="-m 4G",
    threads: 8
    wrapper:
        get_wrapper("samtools", "sort")


rule filterChrM:
    input:
        rules.shift_sort.output,
    output:
        opj(OUTDIR, "dedup/{sample}/{sample}.filtered.bam"),
    conda:
        lambda wildcards: get_environment("samtools", "view")
    shell:
        "(samtools view -h {input} | grep -v chrM | samtools view -bS -F 0x4 -q 30 | samtools sort -@ 4 -o {output}) "
        "&& (samtools index {output})"

if DOWNSAMPLE:
    rule downsample:
        input:
            rules.filterChrM.output
        output:
            opj(OUTDIR, "dedup/{sample}/{sample}.sampled.bam"),
        conda:
            lambda wildcards: get_environment("samtools", "view")
        params:
            counts = DOWNSAMPLE
        threads: 8
        shell:
            """
            total=$(samtools view -c {input})
            frac=$(echo "scale=4; {params.counts} / $total" | bc)
            samtools view -bs $frac {input} | samtools sort -@ {threads} -o {output}
            samtools index {output}
            """