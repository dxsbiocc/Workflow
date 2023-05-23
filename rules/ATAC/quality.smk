rule bamCoverage:
    input:
        # Required input, must index
        'dedup/{sample}/{sample}.filtered.bam',
    output:
        # Required output.
        # Output file format should be one of ['bw', 'bigwig', 'bigWig', 'bedgraph', 'bedGraph'].
        'macs2/bigwig/{sample}.norm.bw'
    params:
        # Optional parameters.
        extra = '--binSize 10 --normalizeUsing RPGC --effectiveGenomeSize ' + str(CHROM_SIZE),
    threads: 1
    log: 
        'logs/deeptools_bamcoverage_{sample}.log'
    wrapper:
        get_wrapper('deeptools', 'bamcoverage')

rule TSSEnrichment:
    input:
        # bigwig = "macs2/bigwig/{pair}.norm.bw", 
        bigwig = get_bigwig,
        bed = config['data']['gtf']
    output:
        matrix_gz = "macs2/matrix/{control}.tss.matrix.gz",   # required
    log:
        "logs/deeptools_compute_tss_matrix_{control}.log"
    params:
        # required argument, choose "scale-regions" or "reference-point"
        command = "reference-point",
        # optional parameters
        extra = "-a 3000 -b 3000 -p 1",
        labels = lambda wildcards, input: [f.split('/')[-1].strip('.norm.bw') for f in input.bigwig]
    wrapper:
        get_wrapper("deeptools", "computematrix")

rule plotTSSHeatmap:
    input:
         rules.TSSEnrichment.output.matrix_gz
    output:
        heatmap_img="macs2/matrix/{control}.tss.heatmap.pdf",  # required
    log:
        "logs/deeptools_TSS_heatmap_{control}.log"
    params:
        # optional parameters
        extra = "--plotType=fill --colorMap Reds Blues ",
    wrapper:
        get_wrapper("deeptools", "plotheatmap")

rule plotTSSProfile:
    input:
        rules.TSSEnrichment.output.matrix_gz
    output:
        plot_img="macs2/matrix/{control}.tss.profile.pdf",  # required
    log:
        "logs/deeptools_TSS_profile_{control}.log"
    params:
        # optional parameters
        extra = "--plotType=fill "
        "--perGroup "
        # "--colors red yellow blue "
        "--dpi 150 ",
    wrapper:
        get_wrapper("deeptools", "plotprofile")

rule genbodyEnrichment:
    input:
        bigwig = get_bigwig,
        bed = config['data']['gtf']
    output:
        matrix_gz = "macs2/matrix/{control}.genebody.matrix.gz",   # required
    log:
        "logs/deeptools_compute_genebody_matrix_{control}.log"
    params:
        # required argument, choose "scale-regions" or "reference-point"
        command="scale-regions",
        # optional parameters
        extra="-a 3000 -b 3000 -p 1 --regionBodyLength 5000 --skipZeros --missingDataAsZero",
        labels = lambda wildcards, input: [f.split('/')[-1].strip('.norm.bw') for f in input.bigwig]
    wrapper:
        get_wrapper("deeptools", "computematrix")

rule plotBodyHeatmap:
    input:
         rules.genbodyEnrichment.output.matrix_gz
    output:
        heatmap_img = "macs2/matrix/{control}.genebody.heatmap.pdf",  # required
    log:
        "logs/deeptools_genebody_heatmap_{control}.log"
    params:
        # optional parameters
        extra = "--plotType=fill --colorMap Reds Blues ",
    wrapper:
        get_wrapper("deeptools", "plotheatmap")

rule ploGenebodytProfile:
    input:
        rules.genbodyEnrichment.output.matrix_gz
    output:
        plot_img="macs2/matrix/{control}.genebody.profile.pdf",  # required
    log:
        "logs/deeptools_genebody_profile_{control}.log"
    params:
        # optional parameters
        extra = "--plotType=fill "
        "--perGroup "
        # "--colors red yellow blue "
        "--dpi 150 ",
    wrapper:
        get_wrapper("deeptools", "plotprofile")



rule multibigwigsummary:
    input:
        # Required input.
        # bigwig = get_bigwig,
        bigwig = expand("macs2/bigwig/{sample}.norm.bw", sample=SAMPLES),
    output:
        # Required output.
        out = 'macs2/bigwig/scores_per_bin.npz',
        # raw_count = ''
    params:
        # Optional parameters.
        subcommand = 'bins',
        extra = '',
        labels = lambda wildcards, input: [f.split('/')[-1].strip('.norm.bw') for f in input.bigwig]
    threads: 1
    log: 
        'logs/deeptools_multibigwigsummary.log'
    wrapper: 
        get_wrapper("deeptools", "multibigwigsummary")
        
rule plotcorrelation:
    input:
        # Required input.
        rules.multibigwigsummary.output.out,
    output:
        img = 'macs2/bigwig/heatmap_spearman_corr_readCounts.pdf',
        # matrix = 'SpearmanCorr_readCounts.tab'
    params:
        # Optional parameters.
        extra = '--skipZeros --whatToPlot heatmap --colorMap RdYlBu --plotNumbers',
        cor = 'spearman',
        title = 'Spearman Correlation of Read Counts',
    threads: 1
    log: 
        'logs/deeptools_plotcorrelation.log'
    wrapper: 
        get_wrapper("deeptools", "plotcorrelation")