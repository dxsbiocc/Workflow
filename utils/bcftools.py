import sys


def infer_out_format(output, uncompressed_bcf=False):
    if output.endswith(".vcf"):
        return "v"
    elif output.endswith(".vcf.gz"):
        return "z"
    elif output.endswith(".bcf"):
        if uncompressed_bcf:
            return "u"
        else:
            return "b"
    else:
        raise ValueError("invalid output file format ('.vcf', '.vcf.gz', '.bcf').")


def get_bcftools_opts(
    snakemake,
    parse_threads=True,
    parse_ref=True,
    parse_region=True,
    parse_samples=True,
    parse_targets=True,
    parse_output=True,
    parse_output_format=True,
    parse_memory=True,
):
    """Obtain bcftools_opts from output, params, and handle resource definitions."""
    bcftools_opts = ""
    extra = snakemake.params.get("extra", "")

    ###############
    ### Threads ###
    ###############
    if parse_threads:
        if "--threads" in extra:
            sys.exit(
                "You have specified number of threads (`--threads`) in `params.extra`; please use `threads`."
            )
        bcftools_opts += (
            ""
            if snakemake.threads <= 1
            else "--threads {}".format(snakemake.threads - 1)
        )

    ######################
    ### Reference file ###
    ######################
    if parse_ref:
        if "-f" in extra or "--fasta-ref" in extra:
            sys.exit(
                "You have specified reference file (`-f/--fasta-ref`) in `params.extra`; this is automatically infered from the `ref` input file."
            )

        if snakemake.input.get("ref"):
            bcftools_opts += f" --fasta-ref {snakemake.input.ref}"

    ####################
    ### Regions file ###
    ####################
    if parse_region:
        if "--region-file" in extra or "-R" in extra:
            sys.exit(
                "You have specified region file (`-R/--regions-file`) in `params.extra`; this is automatically infered from the `regions` input file."
            )

        if snakemake.input.get("regions"):
            bcftools_opts += f" --regions-file {snakemake.input.regions}"

    ####################
    ### Samples file ###
    ####################
    if parse_samples:
        if "-S" in extra or "--samples-file" in extra:
            sys.exit(
                "You have specified samples file (`-S/--samples-file`) in `params.extra`; this is automatically infered from the `samples` input file."
            )

        if snakemake.input.get("samples"):
            bcftools_opts += f" --samples-file {snakemake.input.samples}"

    ####################
    ### Targets file ###
    ####################
    if parse_targets:
        if "-T" in extra or "--targets-file" in extra:
            sys.exit(
                "You have specified samples file (`-T/--targets-file`) in `params.extra`; this is automatically infered from the `targets` input file."
            )

        if snakemake.input.get("targets"):
            bcftools_opts += f" --targets-file {snakemake.input.targets}"

    ###################
    ### Output file ###
    ###################
    if parse_output:
        if "-o" in extra or "--output" in extra:
            sys.exit(
                "You have specified output file (`-o/--output`) in `params.extra`; this is automatically infered from the first output file."
            )
        bcftools_opts += f" -o {snakemake.output[0]}"

    #####################
    ### Output format ###
    #####################
    if parse_output_format:
        if "-O" in extra or "--output-type" in extra:
            sys.exit(
                "You have specified output format (`-O/--output-type`) in `params.extra`; this is automatically infered from output file extension."
            )

        out_format = infer_out_format(
            snakemake.output[0], snakemake.params.get("uncompressed_bcf", False)
        )
        bcftools_opts += f" --output-type {out_format}"

    ##############
    ### Memory ###
    ##############
    if parse_memory:
        if "-m" in extra or "--max-mem" in extra:
            sys.exit(
                "You have provided `-m/--max-mem` in `params.extra`; please use `resources.mem_mb`."
            )
        # Getting memory in megabytes, as advised in documentation.
        if "mem_mb" in snakemake.resources.keys():
            bcftools_opts += " --max-mem {}M".format(snakemake.resources["mem_mb"])
        # Getting memory in gigabytes, for user convenience. Please prefer the use
        # of mem_mb over mem_gb as advised in documentation.
        elif "mem_gb" in snakemake.resources.keys():
            bcftools_opts += " --max-mem {}G".format(snakemake.resources["mem_gb"])

    ################
    ### Temp dir ###
    ################
    if "-T" in extra or "--temp-dir" in extra or "--temp-prefix" in extra:
        sys.exit(
            "You have provided `-T/--temp-dir/--temp-prefix` in `params.extra`; please use `resources.tmpdir`."
        )

    return bcftools_opts
