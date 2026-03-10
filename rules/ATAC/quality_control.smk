# compute the reads count distribution of chromatin
rule plot_chromosome_coverage:
    input:
        rules.dedup_idxstats.output
    output:
        opj(OUTDIR, "QC/chromosome_coverage/{sample}.pdf")
    log:
        opj(OUTDIR, "logs/qc/plot_chromosome_coverage_{sample}.log")
    params:
        width = 10,
        height = 6
    wrapper:
        get_wrapper('scripts', 'R', 'plot_chromosome_coverage')

# frip
rule frip:
    input:
        bam = opj(OUTDIR, "dedup/{pair}/{pair}.filtered.bam") if not DOWNSAMPLE \
         else opj(OUTDIR, "dedup/{pair}/{pair}.sampled.bam"),
        bed = get_peakfile
    output:
        opj(OUTDIR, "QC/frip/{pair}.txt")
    log:
        opj(OUTDIR, "logs/qc/frip_{pair}.log")
    wrapper:
        get_wrapper('scripts', 'Python', 'frip')

# preseq
rule preseq:
    input:
        opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam"),
    output:
        opj(OUTDIR, "QC/preseq/{sample}.lc_extrap")
    log:
        opj(OUTDIR, "logs/qc/preseq_{sample}.log")
    params:
        extra = "-v -seed 1234",
        paired = PAIRED
    wrapper:
        get_wrapper('preseq', 'lc_extrap')

rule plot_preseq:
    input:
        expand(opj(OUTDIR, "QC/preseq/{sample}.lc_extrap"), sample=SAMPLES)
    output:
        opj(OUTDIR, "QC/preseq/{sample}.pdf")
    log:
        opj(OUTDIR, "logs/qc/preseq_plot_{sample}.log")
    params:
        width = 5,
        height = 4,
    wrapper:
        get_wrapper('scripts', 'R', 'plot_preseq')

# fraction of reads in annotated regions
rule frac_anno_regions:
    input:
        bam = opj(OUTDIR, "dedup/{sample}/{sample}.filtered.bam")
    output:
        opj(OUTDIR, "QC/anno_regions/{sample}.txt")
    log:
        opj(OUTDIR, "logs/qc/anno_regions_{sample}.log")
    params:
        dnase = config['parameters']['anno_regions']['dnase'],
        blacklist = config['parameters']['anno_regions']['blacklist'],
        promoter = config['parameters']['anno_regions']['promoter'],
        enhancer = config['parameters']['anno_regions']['enhancer'],
    wrapper:
        get_wrapper('scripts', 'Python', 'anno_regions')

# fragment length stat
rule collectinsertsizemetrics:
    input:
        rules.markduplicates.output.bam
    output:
        txt = opj(OUTDIR, "QC/insertsize/{sample}.hist_data.txt"),
        pdf = opj(OUTDIR, "QC/insertsize/{sample}.hist_graph.pdf"),
    log:
        opj(OUTDIR, "logs/picard/collectinsertsizemetrics_{sample}.log"),
    params:
        extra = "--VERBOSITY ERROR --QUIET true --USE_JDK_DEFLATER true --USE_JDK_INFLATER true --HISTOGRAM_WIDTH 1000 --STOP_AFTER 5000000",
    resources:
        mem_mb = 4096,
    wrapper:
        get_wrapper('picard', 'collectinsertsizemetrics')

rule plot_insertsizemetrics:
    input:
        rules.collectinsertsizemetrics.output.txt
    output:
        qc = opj(OUTDIR, "QC/insertsize/{sample}.isize_dist.txt"),
        pdf = opj(OUTDIR, "QC/insertsize/{sample}.isize_dist.pdf")
    log:
        opj(OUTDIR, "logs/qc/insertsizemetrics_{sample}.log")
    wrapper:
        get_wrapper('scripts', 'Python', 'plot_insertsizemetrics')

# GC bias
rule collectgcbiasmetrics:
    input:
        # BAM aligned to reference genome
        bam = opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam"),
        # reference genome FASTA from which GC-context is inferred
        ref = REFERENCE
    output:
        metrics = opj(OUTDIR, "QC/gcbias/{sample}.gcmetrics.txt"),
        chart = opj(OUTDIR, "QC/gcbias/{sample}.plot.pdf"),
        summary = opj(OUTDIR, "QC/gcbias/{sample}.summary.txt"),
    log:
        opj(OUTDIR, "logs/picard/collectgcbiasmetrics_{sample}.log"),
    params:
        extra = "--USE_JDK_DEFLATER true --USE_JDK_INFLATER true --VERBOSITY ERROR --QUIET true --ASSUME_SORTED false",
    resources:
        mem_mb = 4096,
    wrapper:
        get_wrapper('picard', 'collectgcbiasmetrics')

rule plot_gcbiasmetrics:
    input:
        rules.collectgcbiasmetrics.output.metrics
    output:
        opj(OUTDIR, "QC/gcbias/{sample}.fraglen_dist.pdf")
    log:
        opj(OUTDIR, "logs/qc/gcbiasmetrics_{sample}.log")
    params:
        width = 10,
        height = 6,
    wrapper:
        get_wrapper('scripts', 'R', 'plot_gcbiasmetrics')

# pbc_qc
rule lib_complexity:
    input:
        opj(OUTDIR, "mapped/{sample}/{sample}.sorted.bam"),
    output:
        opj(OUTDIR, "QC/lib_complexity/{sample}.txt")
    log:
        opj(OUTDIR, "logs/qc/lib_complexity_{sample}.log")
    params:
        paired = PAIRED,
        mito_name = 'chrM'
    wrapper:
        get_wrapper('scripts', 'Python', 'lib_complexity')

rule qc_summary:
    input:
        mapping_files = expand(opj(OUTDIR, "logs/mapped/bowtie2_{sample}{suffix}.summary"), sample=SAMPLES, suffix=['', '_spikein']),
        dedup_fils = expand(opj(OUTDIR, "dedup/{sample}/{sample}.metrics.txt"), sample=SAMPLES),
        frip_files = expand(opj(OUTDIR, "QC/frip/{pair}.txt"), pair=PAIRS),
        isize_files = expand(opj(OUTDIR, "QC/anno_regions/{sample}.txt"), sample=SAMPLES),
        lib_complexity_files = expand(opj(OUTDIR, "QC/lib_complexity/{sample}.txt"), sample=SAMPLES),
        idxstats_files = expand(opj(OUTDIR, "QC/mapped/stats/{sample}.idxstats"), sample=SAMPLES),
        anno_files = expand(opj(OUTDIR, "macs2/anno/{pair}.peakAnno.txt"), pair=PAIRS)
    output:
        opj(OUTDIR, "QC/QC_Summary.xlsx")
    log:
        opj(OUTDIR, "logs/qc/summary.log")
    params:
        paired = PAIRED,
        mapping = MAPPING,
        dedup = DEDUP,
    wrapper:
        get_wrapper('scripts', 'Python', 'qc_summary')