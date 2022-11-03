rule shift:
    input:
         rules.mark_duplicates.output.bam,
         rules.cons_index.output
    output:
        "dedup/{sample}.shift.bam"
    log:
        "logs/dedup_shift_{sample}.log"
    threads:
        4
    params:
        # optional parameters
        extra = "{shift} --samFlagInclude {flag} --blackListFileName {bl} --minMappingQuality {mapq}".format(
            shift=config['parameters']['alignmentSieve']['shift'],
            flag=config['parameters']['alignmentSieve']['flag'],
            bl=BLACKLIST,
            mapq=config['parameters']['alignmentSieve']['mapq'],
        )
    wrapper:
        get_wrapper("deeptools", "alignmentSieve")

rule shift_sort:
    input:
        rules.shift.output
    output:
        "dedup/{sample}.shift.sort.bam"
    log:
        "logs/dedup_shift_sort_{sample}.log",
    params:
        extra="-m 4G",
    threads: 8
    wrapper:
        get_wrapper('samtools', 'sort')

rule filterChrM:
    input:
        rules.shift_sort.output,
    output:
        "dedup/{sample}.filtered.bam"
    run:
        shell("(samtools view -h {input} | grep -v chrM | samtools view -bS | samtools sort -@ 4 -o {output})")
        shell("(samtools index {output})")
