rule humann:
    input:
        fastq = opj(OUTDIR, "downsample/{sample}.sampled.fq.gz"),
    output:
        directory(opj(OUTDIR, "humann/{sample}"))
    log:
        opj(OUTDIR, "logs/humann/{sample}.humann.log"),
    params:
        extra = config['parameters']['humann'].get('extra', ''),
        nucleotide = config['parameters']['humann']['nucleotide'],
        protein = config['parameters']['humann']['protein'],
        utility_mapping = config['parameters']['humann']['utility_mapping'],
        metaphlan = f"--bowtie2db {config['parameters']['humann']['metaphlan'].get('db')} "
                    f"--index {config['parameters']['humann']['metaphlan'].get('index')} "
                    f"{config['parameters']['humann']['metaphlan'].get('extra')}"
    threads: config['parameters']['humann']['threads']
    wrapper:
        get_wrapper("humann", "basic")

rule humann_taxonomic:
    input:
        expand(opj(OUTDIR, "humann/{sample}"), sample=SAMPLES)
    output:
        tsv = opj(OUTDIR, "humann/metaphlan4/taxonomy.tsv"),
        spf = opj(OUTDIR, "humann/metaphlan4/taxonomy.classified.spf"),
    log:
        opj(OUTDIR, "logs/humann/metaphlan4/taxonomic.log")
    wrapper:
        get_wrapper("humann", "taxonomic")

rule humann_join:
    input:
        opj(OUTDIR, "humann/")
    output:
        opj(OUTDIR, "humann/abundance.tsv"),
    log:
        opj(OUTDIR, "logs/humann/join_abundance.log")
    params:
        extra = "--file_name pathabundance -s",
    wrapper:
        get_wrapper("humann", "join")

rule humann_renorm_abundance:
    input:
        rules.humann_join.output,
    output:
        opj(OUTDIR, "humann/abundance.relative.tsv"),
    log:
        opj(OUTDIR, "logs/humann/renorm_abundance.log")
    params:
        extra = "--units relab",
    wrapper:
        get_wrapper("humann", "normalization")

rule humann_split_abundance:
    input:
        rules.humann_renorm_abundance.output,
    output:
        stratified = opj(OUTDIR, "humann/abundance.stratified.tsv"),
        unstratified = opj(OUTDIR, "humann/abundance.unstratified.tsv"),
    log:
        opj(OUTDIR, "logs/humann/split_abundance.log")
    wrapper:
        get_wrapper("humann", "split")

rule humann_associate:
    input:
        abundance = opj(OUTDIR, "humann/abundance.tsv"),
    output:
        pcl = opj(OUTDIR, "humann/abundance.pcl"),
        associate = opj(OUTDIR, "humann/associate.tsv"),
    log:
        opj(OUTDIR, "logs/humann/associate.log")
    params:
        group = GROUPS_LINE
    wrapper:
        get_wrapper("humann", "associate")

rule humann_barplot:
    input:
        pcl = opj(OUTDIR, "humann/abundance.pcl"),
    output:
        directory(opj(OUTDIR, "humann/barplot")),
    log:
        opj(OUTDIR, "logs/humann/barplot.log")
    params:
        extra = "--focal-metadata Group --last-metadata Group --sort sum metadata",
        top = config['parameters']['humann'].get('barplot', 10),
    wrapper:
        get_wrapper("humann", "barplot")

rule humann_regroup_ko:
    input:
        opj(OUTDIR, "humann/{sample}"),
    output:
        ko = opj(OUTDIR, "humann/ko/{sample}.ko.tsv"),
    log:
        opj(OUTDIR, "logs/humann/{sample}/regroup_ko.log")
    params:
        extra = "-g uniref90_ko",
    wrapper:
        get_wrapper("humann", "regroup")

rule humann_join_ko:
    input:
        expand(opj(OUTDIR, "humann/ko/{sample}.ko.tsv"), sample=SAMPLES)
    output:
        opj(OUTDIR, "humann/ko/ko.tsv"),
    log:
        opj(OUTDIR, "logs/humann/join_ko.log")
    params:
        extra = "--file_name ko",
    wrapper:
        get_wrapper("humann", "join")


rule humann_split_ko:
    input:
        rules.humann_join_ko.output,
    output:
        stratified = opj(OUTDIR, "humann/ko/ko.stratified.tsv"),
        unstratified = opj(OUTDIR, "humann/ko/ko.unstratified.tsv"),
    log:
        opj(OUTDIR, "logs/humann/split_ko.log")
    wrapper:
        get_wrapper("humann", "split")

rule graphlan_export:
    input:
        taxonomy = rules.humann_taxonomic.output.tsv,
    output:
        tree = opj(OUTDIR, "graphlan/tree.txt"),
        annotation = opj(OUTDIR, "graphlan/annotation.txt"),
    log:
        opj(OUTDIR, "logs/graphlan/export.log")
    params:
        extra = config['parameters']['graphlan']['extra'],
    wrapper:
        get_wrapper("graphlan", "export")

rule graphlan_annotate:
    input:
        tree = rules.graphlan_export.output.tree,
        annotation = rules.graphlan_export.output.annotation,
    output:
        annotated_tree = opj(OUTDIR, "graphlan/annotated_tree.xml"),
    log:
        opj(OUTDIR, "logs/graphlan/annotate.log")
    wrapper:
        get_wrapper("graphlan", "annotate")

rule graphlan_graphlan:
    input:
        rules.graphlan_annotate.output.annotated_tree,
    output:
        opj(OUTDIR, "graphlan/graphlan.pdf"),
    log:
        opj(OUTDIR, "logs/graphlan/graphlan.log")
    params:
        extra = "--external_legends",
    wrapper:
        get_wrapper("graphlan", "graphlan")

rule lefse_format:
    input:
        rules.humann_taxonomic.output.tsv,
    output:
        lefse_input = opj(OUTDIR, "lefse/lefse.txt"),
        lefse_format = opj(OUTDIR, "lefse/input.in"),
    log:
        opj(OUTDIR, "logs/lefse/format.log")
    params:
        extra = "-c 1 -s 2 -u 3 -o 1000000",
        group_line = f"Group\t{GROUPS_LINE}\n"
    wrapper:
        get_wrapper("lefse", "format")

rule lefse_run:
    input:
        rules.lefse_format.output.lefse_format,
    output:
        opj(OUTDIR, "lefse/result.txt"),
    log:
        opj(OUTDIR, "logs/lefse/run.log")
    params:
        extra = "",
    wrapper:
        get_wrapper("lefse", "run")

rule lefse_plot_cladogram:
    input:
        rules.lefse_run.output
    output:
        opj(OUTDIR, "lefse/cladogram.pdf"),
    log:
        opj(OUTDIR, "logs/lefse/cladogram.log")
    params:
        extra = "--format pdf",
    wrapper:
        get_wrapper("lefse", "cladogram")

rule lefse_plot_res:
    input:
        rules.lefse_run.output
    output:
        opj(OUTDIR, "lefse/result.pdf"),
    log:
        opj(OUTDIR, "logs/lefse/result.log")
    params:
        extra = "--format pdf",
    wrapper:
        get_wrapper("lefse", "res")

rule lefse_plot_features:
    input:
        dataset = rules.lefse_format.output.lefse_format,
        lefse = rules.lefse_run.output
    output:
        directory(opj(OUTDIR, "lefse/features")),
    log:
        opj(OUTDIR, "logs/lefse/features.log")
    params:
        extra = "-f diff --archive none --format pdf",
    wrapper:
        get_wrapper("lefse", "features")

rule kraken2:
    input:
        fastq = expand(opj(OUTDIR, "trimmed/{{sample}}/{{sample}}.clean.{run}.fq.gz"), run=RUN),
        db = config['parameters']['kraken2']['db'],
    output:
        output = opj(OUTDIR, "kraken2/{sample}/{sample}.output"),
        report = opj(OUTDIR, "kraken2/{sample}/{sample}.report"),
    log:
        opj(OUTDIR, "logs/kraken2/{sample}.log")
    params:
        extra = "--use-names --report-zero-counts",
    wrapper:
        get_wrapper("kraken2")

rule krakentools_kreport2mpa:
    input:
        rules.kraken2.output.report,
    output:
        opj(OUTDIR, "kraken2/{sample}/{sample}.mpa"),
    log:
        opj(OUTDIR, "logs/kraken2/{sample}.kreport2mpa.log")
    params:
        extra = "--display-header",
        subcommand = "kreport2mpa",
    wrapper:
        get_wrapper("krakentools", "kreport")

rule krakentools_kreport2mpa_percentages:
    input:
        rules.kraken2.output.report,
    output:
        opj(OUTDIR, "kraken2/{sample}/{sample}.percentages.mpa"),
    log:
        opj(OUTDIR, "logs/kraken2/{sample}.kreport2mpa.percentages.log")
    params:
        extra = "--display-header --percentages ",
        subcommand = "kreport2mpa",
    wrapper:
        get_wrapper("krakentools", "kreport")

rule merge_mpa:
    input:
        expand(opj(OUTDIR, "kraken2/{sample}/{sample}.mpa"), sample=SAMPLES),
    output:
        opj(OUTDIR, "kraken2/merged.mpa.tsv"),
    log:
        opj(OUTDIR, "logs/kraken2/merge_mpa.log")
    params:
        usecols = [0, 1],
        column_name = SAMPLES,
        header = 'infer',
    wrapper:
        get_wrapper("scripts", "Python", "merge_samples")

rule bracken:
    input:
        kraken2_report = rules.kraken2.output.report,
        db = config['parameters']['kraken2']['db'],
    output:
        output = opj(OUTDIR, "bracken/{taxonomy}/{sample}.bracken"),
        bracken_report = opj(OUTDIR, "bracken/{taxonomy}/{sample}.report"),
    wildcard_constraints:
        taxonomy = "|".join(TAXONOMY),
    log:
        opj(OUTDIR, "logs/bracken/{taxonomy}/{sample}.log")
    params:
        extra = f"-t {config['parameters']['bracken']['threshold']} " \
                f"-r {config['parameters']['bracken']['readLength']} ",
    wrapper:
        get_wrapper("bracken")

rule merge_bracken:
    input:
        expand(opj(OUTDIR, "bracken/{{taxonomy}}/{sample}.bracken"), sample=SAMPLES),
    output:
        opj(OUTDIR, "bracken/{taxonomy}/merged.bracken.tsv"),
    log:
        opj(OUTDIR, "logs/bracken/{taxonomy}/merge_bracken.log")
    params:
        usecols = [0, 5],
        header = 'infer',
        column_name = SAMPLES,
    wildcard_constraints:
        taxonomy = "|".join(TAXONOMY),
    wrapper:
        get_wrapper("scripts", "Python", "merge_samples")

rule krakentools_alpha_diversity:
    input:
        expand(opj(OUTDIR, "bracken/{{taxonomy}}/{sample}.bracken"), sample=SAMPLES),
    output:
        opj(OUTDIR, "bracken/{taxonomy}/alpha.tsv"),
    log:
        opj(OUTDIR, "logs/bracken/{taxonomy}/alpha_diversity.log")
    params:
        metrics = ["Sh", "BP", "Si", "ISi", "F"],
        sample_list = SAMPLES,
    wildcard_constraints:
        taxonomy = "|".join(TAXONOMY),
    wrapper:
        get_wrapper("krakentools", "alpha")

rule krakentools_beta_diversity:
    input:
        expand(opj(OUTDIR, "bracken/{{taxonomy}}/{sample}.bracken"), sample=SAMPLES),
    output:
        opj(OUTDIR, "bracken/{taxonomy}/beta.txt"),
    log:
        opj(OUTDIR, "logs/bracken/{taxonomy}/beta_diversity.log")
    params:
        type = "bracken",
        level = "all",
        extra = "",
    wildcard_constraints:
        taxonomy = "|".join(TAXONOMY),
    wrapper:
        get_wrapper("krakentools", "beta")

rule krakentools_kreport2krona:
    input:
        rules.bracken.output.bracken_report,
    output:
        opj(OUTDIR, "bracken/{taxonomy}/{sample}.krona"),
    log:
        opj(OUTDIR, "logs/bracken/{taxonomy}/{sample}.kreport2krona.log")
    params:
        extra = "--no-intermediate-ranks",
        subcommand = "kreport2krona",
    wrapper:
        get_wrapper("krakentools", "kreport")

rule krakentools_ktImportText:
    input:
        rules.krakentools_kreport2krona.output,
    output:
        opj(OUTDIR, "bracken/{taxonomy}/{sample}.krona.html"),
    log:
        opj(OUTDIR, "logs/bracken/{taxonomy}/{sample}.ktImportText.log")
    params:
        extra = "",
    wrapper:
        get_wrapper("krona", "importtext")
