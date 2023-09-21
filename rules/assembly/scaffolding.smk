rule download_juicer_tools:
    output:
        paths.juicebox.tools,
    params:
        url = config["software"]["juicer"]["tools_url"],
    conda:
        "../envs/juicer_tools.yml"
    log:
        to_log(paths.juicebox.tools),
    shell:
        "wget -O - {params.url} > {output} 2>{log}"

rule bwa_mem2_index:
    input:
        "genome/genome.{hap}.fasta",
    output:
        multiext("scaffolding/index/{hap}", ".0123", ".amb", ".ann", ".pac", ".bwt.2bit.64"),
    log:
        "logs/bwa-mem2/index.log",
    threads: 20
    wrapper:
        get_wrapper("bwa-mem2", "index")

rule hicfindrestsite:
    input:
        "genome/genome.{hap}.fasta"
    output:
        "scaffolding/restsite/DpnII.bed"
    log:
        "logs/hicexplorer/hicfindrestsite.log"
    params:
        pattern = ["GATC"],  # HindIII: 'AAGCTT', DpnII/MboI/Sau3AI: 'GATC', Arima: 'GATC', 'GANTC'
        extra = ""
    threads: 20
    wrapper:
        get_wrapper("hicexplorer", "pre-processing", "hicfindrestsite")

rule bwa_mem2_mem:
    input:
        reads = expand("trimmed/{{hic}}/{{hic}}.clean.{run}.fq.gz", run=RUN),
        # Index can be a list of (all) files created by bwa, or one of them
        idx = multiext("scaffolding/index/{hap}", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    output:
        "scaffolding/mapped/{hic}.bam",
    log:
        "logs/bwa_mem2/{hic}.log",
    params:
        extra = r"-R '@RG\tID:{hic}\tSM:{hic}'",
        sort = "samtools",  # Can be 'none', 'samtools', or 'picard'.
        sort_order = "coordinate",  # Can be 'coordinate' (default) or 'queryname'.
        sort_extra = "",  # Extra args for samtools/picard sorts.
    threads: 20
    wrapper:
        get_wrapper("bwa-mem2", "mem")


rule pretext_map:
    input:
        rules.bwa_mem2_mem.output,
    output:
        "scaffolding/pretext/{hic}/map.pretext",
    log:
        "logs/scaffolding/pretext_map_{hic}.log",
    params:
        extra = "--sortby length --sortorder descend --mapq 10",
    wrapper:
        get_wrapper("pretext", "map")


rule pretext_snapshot_png:
    input:
        rules.pretext_map.output,
    output:
        all = directory("scaffolding/pretext/{hic}/all_map"),
        full = "scaffolding/pretext/{hic}/full_map.png",
    log:
        "logs/scaffolding/pretext_snapshot_png_{hic}.log",
    params:
        extra="--resolution 1080",
    wrapper:
        get_wrapper("pretext", "snapshot")