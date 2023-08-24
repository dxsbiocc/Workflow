rule bam2fastq:
    input:
        files = config["reads"].split(" ")
    output:
        "{sample}_combined_reads.fastq.gz"
    params:
        prefix = "{sample}_combined_reads"
    wrapper:
        get_wrapper("bedtools", "bam2fastq")


rule adapterfilt:
    input:
        "{sample}_combined_reads.fastq.gz"
    output:
        temp(hifi_outdir + "/{sample}.filt.fastq.gz")
    conda:
        "envs/adapter.yaml"
    shell:
        "cutadapt --discard-trimmed --overlap=35 -b ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT -b ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT -o {output} {input}"
