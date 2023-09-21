rule find_valid_interactions:
    input:
        "results/aligned_reads.bam",
        fragments
    output:
        "results/valid_interactions.txt"
    shell:
        "hicFindTADs --matrix-format bam --threads 4 --restriction {fragments} {input} {output}"

rule create_interactions_matrix:
    input:
        "results/valid_interactions.txt"
    output:
        "results/interactions_matrix.txt"
    shell:
        "hicConvertFormat --outFileName {output} --inputFormat interactions --coord1=BP1 --coord2=BP2 --value=COUNTS --bedFile={input}"

rule normalize_matrix:
    input:
        "results/interactions_matrix.txt"
    output:
        "results/normalized_matrix.txt"
    shell:
        "hicNormalize --outFileName {output} --matFormat --normMethod {config[normalization_method]} --inFile {input}"

rule plot_heatmap:
    input:
        "results/normalized_matrix.txt"
    output:
        "results/heatmap.png"
    shell:
        "hicPlotMatrix --colorMap {config[color_map]} --log1p --outFileName {output} --matrix {input}"
