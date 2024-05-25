rule macs2Narrow:
    input:
        unpack(get_paired_bam)
    output:
        # all output-files must share the same basename and only differ by it's extension
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
        multiext(opj(OUTDIR, "macs2/narrow/{pair}"),
                "_peaks.xls",   ### required
                "_treat_pileup.bdg",
                "_control_lambda.bdg",
                ### optional output files
                "_peaks.narrowPeak",
                "_summits.bed"
                )
    log:
        opj(OUTDIR, "logs/macs2/macs2_callpeak_narrow_{pair}.log")
    params:
        extra = lambda wildcards: get_macs2(wildcards.pair, True),
    wrapper:
        get_wrapper('macs2', 'callpeak')

rule macs2Broad:
    input:
        unpack(get_paired_bam)
    output:
        multiext(opj(OUTDIR, "macs2/broad/{pair}"),
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
        opj(OUTDIR, "logs/macs2/macs2_callpeak_broad_{pair}.log")
    params:
        extra = lambda wildcards: get_macs2(wildcards.pair, False)
    wrapper:
        get_wrapper('macs2', 'callpeak')

rule annotatePeaks:
    input:
        peakfile = opj(OUTDIR, "macs2/broad/{pair}_peaks.broadPeak") if PEAKMODE == 'broad' \
                                    else opj(OUTDIR, "macs2/narrow/{pair}_peaks.narrowPeak"),
        gtf = config['data']['gtf']
    output:
        anno_txt = opj(OUTDIR, 'macs2/anno/{pair}.peakAnno.txt'),
        anno_pdf = opj(OUTDIR, 'macs2/anno/{pair}.peakAnno.pdf'),
    log:
        opj(OUTDIR, "logs/anno/peak_anno_{pair}.log")
    wrapper:
        get_wrapper('scripts', 'R', 'annoPeaks')

rule peak2saf:
    input:
        peakfile = opj(OUTDIR, "macs2/broad/{pair}_peaks.broadPeak") if PEAKMODE == 'broad' else \
                   opj(OUTDIR, "macs2/narrow/{pair}_peaks.narrowPeak"),
    output:
        saf = opj(OUTDIR, 'macs2/quantify/{pair}.saf')
    log:
        opj(OUTDIR, "logs/anno/peak2saf_{pair}.log")
    shell:
        """
        awk 'OFS="\t" {{print $4, $1, $2, $3, "."}}' {input.peakfile} > {output.saf}
        """

rule quantifyPeaks:
    input:
        samples = opj(OUTDIR, 'dedup/{pair}/{pair}.filtered.bam'),
        annotation = opj(OUTDIR, 'macs2/quantify/{pair}.saf'),
    output:
        quant = opj(OUTDIR, 'macs2/quantify/{pair}.counts')
    log:
        opj(OUTDIR, "logs/anno/quantifyPeaks_{pair}.log")
    threads: 10
    params:
        extra = '-F SAF',
        paired = PAIRED
    wrapper:
        get_wrapper('subread', 'featurecounts')

rule merge_peaks:
    input:
        expand(opj(OUTDIR, "macs2/broad/{pair}_peaks.broadPeak"), pair=PAIRS) if PEAKMODE == 'broad' else \
        expand(opj(OUTDIR, "macs2/narrow/{pair}_peaks.narrowPeak"), pair=PAIRS)
    output:
        bed = opj(OUTDIR, 'macs2/quantify/merged_peaks.bed'),
        saf = opj(OUTDIR, 'macs2/quantify/merged.saf')
    log:
        opj(OUTDIR, 'logs/anno/merge_peaks.log')
    conda:
        get_environment('bedtools', 'merge')
    shell:
        """
        cat {input} | bedtools sort -i - | bedtools merge -i - > {output.bed}
        awk 'OFS="\t" {print $1":"$2"-"$3, $1, $2, $3, "."}' {output.bed} > {output.saf}
        """

rule quantifyMerged:
    input:
        samples = expand(opj(OUTDIR, 'dedup/{pair}/{pair}.filtered.bam'), pair=PAIRS),
        annotation = opj(OUTDIR, 'macs2/quantify/merged.saf'),
    output:
        quant = opj(OUTDIR, 'macs2/quantify/merged.quants')
    log:
        opj(OUTDIR, "logs/anno/quantifyMerged.log")
    threads: 10
    params:
        extra = '-F SAF',
        paired = PAIRED
    wrapper:
        get_wrapper('subread', 'featurecounts')