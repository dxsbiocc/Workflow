rule macs2Narrow:
    input:
        unpack(get_paired_bam)
    output:
        # all output-files must share the same basename and only differ by it's extension
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
        multiext("macs2/narrow/{pair}",
                 "_peaks.xls",   ### required
                 "_treat_pileup.bdg",
                 "_control_lambda.bdg",
                 ### optional output files
                 "_peaks.narrowPeak",
                 "_summits.bed"
                 )
    log:
        "logs/macs2_callpeak_narrow_{pair}.log"
    params:
        extra = lambda wildcards: get_macs2(wildcards.pair, True),
    wrapper:
        get_wrapper('macs2', 'callpeak')

rule macs2Broad:
    input:
        unpack(get_paired_bam)
    output:
        # all output-files must share the same basename and only differ by it's extension
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
        multiext("macs2/broad/{pair}",
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
        "logs/macs2_callpeak_broad_{pair}.log"
    params:
        extra = lambda wildcards: get_macs2(wildcards.pair, False)
    wrapper:
        get_wrapper('macs2', 'callpeak')

rule annotatePeaks:
    input:
        peakfile = "macs2/narrow/{pair}_peaks.narrowPeak",
        gtf = config['data']['gtf']
    output:
        expand('macs2/anno/{{pair}}.peakAnno.{ext}', ext=['pdf', 'txt'])
    log:
        "logs/peak_anno_{pair}.log"
    params:
        sample = lambda wildcards: '{}'.format(wildcards.pair),
        output = config['workdir'] + "/macs2/anno/"
    wrapper:
        get_wrapper('scripts', 'annoPeaks')