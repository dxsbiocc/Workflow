rule megahit:
    input:
        r1 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R1.fq.gz"),
        r2 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R2.fq.gz"),
    output:
        contigs = opj(OUTDIR, "megahit/{sample}/final.contigs.fasta"),
        json = opj(OUTDIR, "megahit/{sample}/options.json"),
        log = opj(OUTDIR, "megahit/{sample}/log"),
    params:
        extra = config['parameters']['megahit']['extra'],
    threads: config['parameters']['megahit']['threads'],
    resources:
        mem_mb = config['parameters']['megahit']['memory'],
    wrapper:
        get_wrapper("megahit")

rule metaspades:
    input:
        r1 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R1.fq.gz"),
        r2 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R2.fq.gz"),
    output:
        contigs = opj(OUTDIR, "metaspades/{sample}/final.contigs.fasta"),
        scaffolds = opj(OUTDIR, "metaspades/{sample}/final.scaffolds.fasta"),
        dir = directory(opj(OUTDIR, "metaspades/{sample}/intermediate_files")),
    log:
        opj(OUTDIR, "logs/metaspades.{sample}.log"),
    params:
        extra = config['parameters']['metaspades']['extra'],
        k = "auto",
        mode = "meta",  # options: meta, isolate, sc, plasmid, metaplasmid, metaviral, rna, rnaviral, bio, corona, sewage
    threads: config['parameters']['metaspades']['threads'],
    resources:
        mem_mb = config['parameters']['metaspades']['memory'],
    wrapper:
        get_wrapper("spades", "metaspades")

rule quast:
    input:
        fasta = opj(OUTDIR, "{assembly}/{sample}/final.contigs.fasta"),
    output:
        report_html=opj(OUTDIR, "quast/{assembly}/{sample}/report.html"),
        report_pdf=opj(OUTDIR, "quast/{assembly}/{sample}/report.pdf"),
        report_tsv=opj(OUTDIR, "quast/{assembly}/{sample}/report.tsv"),
        treport_tex=opj(OUTDIR, "quast/{assembly}/{sample}/treport.tex"),
        treport_txt=opj(OUTDIR, "quast/{assembly}/{sample}/treport.txt"),
        treport_tsv=opj(OUTDIR, "quast/{assembly}/{sample}/treport.tsv"),
        stats_cum=opj(OUTDIR, "quast/{assembly}/{sample}/stats/cumulative.pdf"),
        stats_gc_plot=opj(OUTDIR, "quast/{assembly}/{sample}/stats/gc.pdf"),
        stats_gc_fasta=opj(OUTDIR, "quast/{assembly}/{sample}/stats/gc_fasta.pdf"),
        stats_nx=opj(OUTDIR, "quast/{assembly}/{sample}/stats/Nx.pdf"),
        icarus=opj(OUTDIR, "quast/{assembly}/{sample}/icarus.html"),
        icarus_viewer=opj(OUTDIR, "quast/{assembly}/{sample}/icarus_viewer.html"),
    log:
        opj(OUTDIR, "logs/quast.{assembly}/{sample}.log"),
    params:
        extra = config['parameters']['quast']['extra'],
    wrapper:
        get_wrapper("quast")

if len(ASSEMBLY) == 2:
    rule quast_both:
        input:
            fasta = expand(opj(OUTDIR, "{assembly}/{sample}/final.contigs.fasta"), assembly=ASSEMBLY, sample=SAMPLES),
        output:
            report_html=opj(OUTDIR, "quast/compare/{sample}/report.html"),
            report_pdf=opj(OUTDIR, "quast/compare/{sample}/report.pdf"),
            report_tsv=opj(OUTDIR, "quast/compare/{sample}/report.tsv"),
            treport_tex=opj(OUTDIR, "quast/compare/{sample}/treport.tex"),
            treport_txt=opj(OUTDIR, "quast/compare/{sample}/treport.txt"),
            treport_tsv=opj(OUTDIR, "quast/compare/{sample}/treport.tsv"),
            stats_cum=opj(OUTDIR, "quast/compare/{sample}/stats/cumulative.pdf"),
            stats_gc_plot=opj(OUTDIR, "quast/compare/{sample}/stats/gc.pdf"),
            stats_gc_fasta=opj(OUTDIR, "quast/compare/{sample}/stats/gc_fasta.pdf"),
            stats_nx=opj(OUTDIR, "quast/compare/{sample}/stats/Nx.pdf"),
            icarus=opj(OUTDIR, "quast/compare/{sample}/icarus.html"),
            icarus_viewer=opj(OUTDIR, "quast/compare/{sample}/icarus_viewer.html"),
        params:
            extra = "--label 'megahit,metaspades'",
        wrapper:
            get_wrapper("quast")

rule metaquast:
    input:
        fasta = opj(OUTDIR, "{assembly}/{sample}/final.contigs.fasta"),
    output:
        report_html=opj(OUTDIR, "metaquast/{assembly}/{sample}/report.html"),
        report_pdf=opj(OUTDIR, "metaquast/{assembly}/{sample}/report.pdf"),
        report_tsv=opj(OUTDIR, "metaquast/{assembly}/{sample}/report.tsv"),
        treport_tex=opj(OUTDIR, "metaquast/{assembly}/{sample}/treport.tex"),
        treport_txt=opj(OUTDIR, "metaquast/{assembly}/{sample}/treport.txt"),
        treport_tsv=opj(OUTDIR, "metaquast/{assembly}/{sample}/treport.tsv"),
        stats_cum=opj(OUTDIR, "metaquast/{assembly}/{sample}/stats/cumulative.pdf"),
        stats_gc_plot=opj(OUTDIR, "metaquast/{assembly}/{sample}/stats/gc.pdf"),
        stats_gc_fasta=opj(OUTDIR, "metaquast/{assembly}/{sample}/stats/gc_fasta.pdf"),
        stats_nx=opj(OUTDIR, "metaquast/{assembly}/{sample}/stats/Nx.pdf"),
        icarus=opj(OUTDIR, "metaquast/{assembly}/{sample}/icarus.html"),
        icarus_viewer=opj(OUTDIR, "metaquast/{assembly}/{sample}/icarus_viewer.html"),
    log:
        opj(OUTDIR, "logs/metaquast/{assembly}/{sample}.log"),
    params:
        extra = config['parameters']['quast']['extra'],
        command = "metaquast",
    wrapper:
        get_wrapper("quast")

rule prodigal:
    input:
        input_file = opj(OUTDIR, "{assembly}/{sample}/final.contigs.fasta"),
    output:
        nuc_file = opj(OUTDIR, "prodigal/{assembly}/{sample}/gene.fa"),
        output_file = opj(OUTDIR, "prodigal/{assembly}/{sample}/gene.gff"),
    log:
        opj(OUTDIR, "logs/prodigal.{assembly}/{sample}.log"),
    params:
        extra = config['parameters']['prodigal']['extra'],
    wrapper:
        get_wrapper("prodigal")

rule concat_all_gene:
    input:
        expand(opj(OUTDIR, "prodigal/{{assembly}}/{sample}/gene.fa"), sample=SAMPLES),
    output:
        opj(OUTDIR, "NR/{assembly}/all_gene.fa"),
    log:
        opj(OUTDIR, "logs/{assembly}/concat_all_gene.log"),
    params:
        command = "concat",
    wrapper:
        get_wrapper("seqkit")

rule complete_gene:
    input:
        rules.prodigal.output.nuc_file,
    output:
        opj(OUTDIR, "prodigal/{assembly}/{sample}/complete_gene.fa"),
    log:
        opj(OUTDIR, "logs/complete_gene/{assembly}/{sample}_complete_gene.log"),
    params:
        command = "grep",
        extra = config['parameters']['seqkit']['complete_gene'],
    wrapper:
        get_wrapper("seqkit")

rule concat_full_gene:
    input:
        expand(opj(OUTDIR, "prodigal/{{assembly}}/{sample}/complete_gene.fa"), sample=SAMPLES),
    output:
        opj(OUTDIR, "NR/{assembly}/gene.fa"),
    log:
        opj(OUTDIR, "logs/{assembly}/concat_full_gene.log"),
    params:
        command = "concat",
    wrapper:
        get_wrapper("seqkit")

rule mmseqs2_easy_cluster:
    input:
        fasta = rules.concat_full_gene.output if FULL_GENE else rules.concat_all_gene.output,
    output:
        cluster = opj(OUTDIR, "NR/{assembly}/nucleotide.fa"),
    log:
        opj(OUTDIR, "logs/mmseqs2_easy_cluster/{assembly}.log"),
    params:
        extra = config['parameters']['mmseqs2']['easy-cluster'],
    threads: config['parameters']['mmseqs2']['threads'],
    wrapper:
        get_wrapper("mmseqs2", "easy-cluster")

rule translate_gene:
    input:
        trim = rules.mmseqs2_easy_cluster.output.cluster,
    output:
        opj(OUTDIR, "NR/{assembly}/protein.fa"),
    log:
        opj(OUTDIR, "logs/seqkit/{assembly}/translate_gene.log"),
    params:
        extra = "",
        command = "translate",
    wrapper:
        get_wrapper("seqkit")

rule protein_no_end:
    input:
        rules.translate_gene.output
    output:
        opj(OUTDIR, "NR/{assembly}/protein_no_end.fa"),
    shell:
        "sed 's/\*//g' {input} > {output}"

rule salmon_index:
    input:
        sequences = rules.mmseqs2_easy_cluster.output.cluster,
    output:
        directory(opj(OUTDIR, "NR/{assembly}/index")),
    log:
        opj(OUTDIR, "logs/salmon_index/{assembly}.log"),
    params:
        extra = config['parameters']['salmon']['index']['extra'],
    threads: config['parameters']['salmon']['index']['threads'],
    wrapper:
        get_wrapper("salmon", "index")

rule salmon_quant:
    input:
        r1 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R1.fq.gz"),
        r2 = opj(OUTDIR, "trimmed/{sample}/{sample}.clean.R2.fq.gz"),
        index = opj(OUTDIR, "NR/{assembly}/index"),
    output:
        quant = opj(OUTDIR, "salmon/{assembly}/{sample}/quant.sf"),
        lib = opj(OUTDIR, "salmon/{assembly}/{sample}/lib_format_counts.json"),
    log:
        opj(OUTDIR, "logs/salmon/{assembly}/{sample}.log"),
    params:
        libtype = "A",
        extra = config['parameters']['salmon']['quant']['extra'],
    threads: config['parameters']['salmon']['quant']['threads'],
    wrapper:
        get_wrapper("salmon", "quant")

rule salmon_merge_tpm:
    input:
        quants = expand(opj(OUTDIR, "salmon/{{assembly}}/{sample}/quant.sf"), sample=SAMPLES),
    output:
        opj(OUTDIR, "salmon/{assembly}/quant.tpm"),
    log:
        opj(OUTDIR, "logs/salmon/{assembly}_merge_tpm.log"),
    params:
        extra = ""
    wrapper:
        get_wrapper("salmon", "quantmerge")

rule salmon_merge_counts:
    input:
        quants = expand(opj(OUTDIR, "salmon/{{assembly}}/{sample}/quant.sf"), sample=SAMPLES),
    output:
        opj(OUTDIR, "salmon/{assembly}/quant.counts"),
    log:
        opj(OUTDIR, "logs/salmon/{assembly}_merge_counts.log"),
    params:
        extra = "--column NumReads"
    wrapper:
        get_wrapper("salmon", "quantmerge")

rule emapper:
    input:
        data_dir = config['parameters']['eggnog']['db'],
        protein = rules.protein_no_end.output,
    output:
        opj(OUTDIR, "eggnog/{assembly}/eggnog.emapper.annotations"),
    log:
        opj(OUTDIR, "logs/eggnog/{assembly}.log"),
    params:
        extra = config['parameters']['eggnog']['extra'],
        prefix = "eggnog",
    wildcard_constraints:
        assembly = "|".join(ASSEMBLY),
    wrapper:
        get_wrapper("eggnog")

rule eggnog_anno_output:
    input:
        rules.emapper.output,
    output:
        opj(OUTDIR, "eggnog/{assembly}/eggnog.annotations"),
    wildcard_constraints:
        assembly = "|".join(ASSEMBLY),
    shell:
        "grep -v '^##' {input} | sed '1 s/^#//' > {output}"

rule rgi_card:
    input:
        fasta = rules.protein_no_end.output,
    output:
        opj(OUTDIR, "rgi/{assembly}/protein.json"),
    params:
        extra = config['parameters']['rgi']['extra'],
        aligner = config['parameters']['rgi']['aligner'],  # BLAST and DIAMOND
        intype = config['parameters']['rgi']['intype'],    # contig, protein and read
    threads: config['parameters']['rgi']['thread']
    wrapper:
        get_wrapper("rgi")

rule rgi_nucleotide:
    input:
        fasta = rules.mmseqs2_easy_cluster.output.cluster,
    output:
        opj(OUTDIR, "rgi/{assembly}/nucleotide.json"),
    params:
        extra = config['parameters']['rgi']['extra'],
        aligner = config['parameters']['rgi']['aligner'],  # BLAST and DIAMOND
        intype = "contig",    # contig, protein and read
    threads: config['parameters']['rgi']['thread']
    wrapper:
        get_wrapper("rgi")