rule bamCoverage:
    input:
        # Required input.
        'dedup/{sample}.filtered.bam',
        # must index
        # 'dedup/{sample}.filtered.bam.bai'
    output:
        # Required output.
        # Output file format should be one of ['bw', 'bigwig', 'bigWig', 'bedgraph', 'bedGraph'].
        'macs2/bigwig/{sample}.norm.bw'
    params:
        # Optional parameters.
        extra = '--binSize 10 --normalizeUsing RPGC --effectiveGenomeSize ' + str(total_chrom_size),
    threads: 1
    log: 
        'logs/{sample}_deeptools_bamcoverage.log'
    wrapper:
        get_wrapper('deeptools', 'bamcoverage')

rule TSSEnrichment:
    input:
        bigwig = expand("macs2/bigwig/{sample}.norm.bw", sample=samples.index),
        bed = config['data']['gtf']
    output:
        matrix_gz = "macs2/matrix/tss.matrix.gz",   # required
    log:
        "logs/deeptools_compute_tss_matrix.log"
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
        heatmap_img="macs2/matrix/tss.heatmap.png",  # required
    log:
        "logs/deeptools_TSS_heatmap.log"
    params:
        # optional parameters
        extra = "--plotType=fill --colorMap Reds "
    wrapper:
        get_wrapper("deeptools", "plotheatmap")

rule plotTSSProfile:
    input:
        rules.TSSEnrichment.output.matrix_gz
    output:
        plot_img="macs2/matrix/tss.profile.png",  # required
    log:
        "logs/deeptools_TSS_profile.log"
    params:
        # optional parameters
        extra = "--plotType=fill "
        "--perGroup "
        "--colors red yellow blue "
        "--dpi 150 "
    wrapper:
        get_wrapper("deeptools", "plotprofile")

rule genbodyEnrichment:
    input:
        bigwig = expand("macs2/bigwig/{sample}.norm.bw", sample=samples.index),
        bed = config['data']['gtf']
    output:
        matrix_gz = "macs2/matrix/genebody.matrix.gz",   # required
    log:
        "logs/deeptools_compute_genebody_matrix.log"
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
        heatmap_img = "macs2/matrix/genebody.heatmap.png",  # required
    log:
        "logs/deeptools_genebody_heatmap.log"
    params:
        # optional parameters
        extra = "--plotType=fill --colorMap Reds "
    wrapper:
        get_wrapper("deeptools", "plotheatmap")

rule ploGenebodytProfile:
    input:
        rules.genbodyEnrichment.output.matrix_gz
    output:
        plot_img="macs2/matrix/genebody.profile.png",  # required
    log:
        "logs/deeptools_genebody_profile.log"
    params:
        # optional parameters
        extra = "--plotType=fill "
        "--perGroup "
        "--colors red yellow blue "
        "--dpi 150 "
    wrapper:
        get_wrapper("deeptools", "plotprofile")