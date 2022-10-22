rule macs2_narrow:
    input:
        treatment = "dedup/{sample}.shift.sort.bam",   # required: treatment sample(s)
    output:
        # all output-files must share the same basename and only differ by it's extension
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
        multiext("macs2/narrow/{sample}",
                 "_peaks.xls",   ### required
                 ### optional output files
                 "_peaks.narrowPeak",
                 "_summits.bed"
                 )
    log:
        "log/macs2_{sample}_callpeak_narrow.log"
    params:
        extra = "-f BAMPE -g hs -q 0.01 -B --SPMR --keep-dup all"
    wrapper:
        get_wrapper('macs2', 'callpeak')

rule macs2_broad:
    input:
        treatment = "dedup/{sample}.shift.sort.bam",   # required: treatment sample(s)
    output:
        # all output-files must share the same basename and only differ by it's extension
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
        multiext("macs2/broad/{sample}",
                 "_peaks.xls",   ### required
                 ### optional output files
                 # these output extensions internally set the --bdg or -B option:
                 "_treat_pileup.bdg",
                 "_control_lambda.bdg",
                 # these output extensions internally set the --broad option:
                 "_peaks.broadPeak",
                 "_peaks.gappedPeak"
                 )
    log:
        "log/macs2_{sample}_callpeak_broad.log"
    params:
        extra = "-f BAMPE -g hs --broad --broad-cutoff 0.1 -B --SPMR --keep-dup all"
    wrapper:
        get_wrapper('macs2', 'callpeak')


rule bamCoverage:
    input:
        # Required input.
        'dedup/{sample}.shift.sort.bam',
    output:
        # Required output.
        # Output file format should be one of ['bw', 'bigwig', 'bigWig', 'bedgraph', 'bedGraph'].
        'macs2/bigwig/{sample}.cpm.norm.bw'
    params:
        # Optional parameters.
        extra = '--binSize 10 --normalizeUsing CPM --effectiveGenomeSize' + '{total_chrom_size}',
    threads: 1
    log: 
        'log/deeptools_bamcoverage_{sample}.log'
    wrapper:
        get_wrapper('deeptools', 'bamcoverage')