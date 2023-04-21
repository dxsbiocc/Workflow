rule bamCoverage:
    input:
        # Required input, must index
        'dedup/{sample}.filtered.bam',
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
        matrix_gz = "macs2/matrix/{pair}.tss.matrix.gz",   # required
    log:
        "logs/deeptools_compute_tss_matrix_{pair}.log"
    params:
        # required argument, choose "scale-regions" or "reference-point"
        command="reference-point",
        # optional parameters
        extra="-a 3000 -b 3000 -p 1"
    wrapper:
        get_wrapper("deeptools", "computematrix")

rule plotTSSHeatmap:
    input:
         rules.TSSEnrichment.output.matrix_gz
    output:
        heatmap_img="macs2/matrix/{pair}.tss.heatmap.pdf",  # required
    log:
        "logs/deeptools_TSS_heatmap_{pair}.log"
    params:
        # optional parameters
        extra = "--plotType=fill --colorMap Reds Blues "
    wrapper:
        get_wrapper("deeptools", "plotheatmap")

rule plotTSSProfile:
    input:
        rules.TSSEnrichment.output.matrix_gz
    output:
        plot_img="macs2/matrix/{pair}.tss.profile.pdf",  # required
    log:
        "logs/deeptools_TSS_profile_{pair}.log"
    params:
        # optional parameters
        extra = "--plotType=fill "
        "--perGroup "
        # "--colors red yellow blue "
        "--dpi 150 "
    wrapper:
        get_wrapper("deeptools", "plotprofile")

rule genbodyEnrichment:
    input:
        bigwig = get_bigwig,
        bed = config['data']['gtf']
    output:
        matrix_gz = "macs2/matrix/{pair}.genebody.matrix.gz",   # required
    log:
        "logs/deeptools_compute_genebody_matrix_{pair}.log"
    params:
        # required argument, choose "scale-regions" or "reference-point"
        command="scale-regions",
        # optional parameters
        extra="-a 3000 -b 3000 -p 1 --regionBodyLength 5000 --skipZeros --missingDataAsZero"
    wrapper:
        get_wrapper("deeptools", "computematrix")

rule plotBodyHeatmap:
    input:
         rules.genbodyEnrichment.output.matrix_gz
    output:
        heatmap_img = "macs2/matrix/{pair}.genebody.heatmap.pdf",  # required
    log:
        "logs/deeptools_genebody_heatmap_{pair}.log"
    params:
        # optional parameters
        extra = "--plotType=fill --colorMap Reds Blues "
    wrapper:
        get_wrapper("deeptools", "plotheatmap")

rule ploGenebodytProfile:
    input:
        rules.genbodyEnrichment.output.matrix_gz
    output:
        plot_img="macs2/matrix/{pair}.genebody.profile.pdf",  # required
    log:
        "logs/deeptools_genebody_profile_{pair}.log"
    params:
        # optional parameters
        extra = "--plotType=fill "
        "--perGroup "
        # "--colors red yellow blue "
        "--dpi 150 "
    wrapper:
        get_wrapper("deeptools", "plotprofile")