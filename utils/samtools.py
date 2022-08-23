import sys
from os import path


def infer_out_format(file_name):
    out_name, out_ext = path.splitext(file_name)
    return out_ext[1:].upper()


def get_samtools_opts(
    snakemake,
    parse_threads=True,
    parse_ref=True,
    parse_write_index=True,
    parse_output=True,
    parse_output_format=True,
):
    """Obtain samtools_opts from output, params, and handle resource definitions in resources."""
    samtools_opts = ""
    extra = snakemake.params.get("extra", "")
    idx = snakemake.output.get("idx", "")

    ###############
    ### Threads ###
    ###############
    if parse_threads:
        if "-@" in extra or "--threads" in extra:
            sys.exit(
                "You have specified number of threads (`-@/--threads`) in params.extra; please use only `threads`."
            )
        samtools_opts += (
            ""
            if snakemake.threads <= 1
            else " --threads {}".format(snakemake.threads - 1)
        )

    ######################
    ### Reference file ###
    ######################
    if parse_ref:
        if "--reference" in extra:
            sys.exit(
                "You have specified reference file (`--reference`) in `params.extra`; this is automatically infered from the `ref` input file."
            )

        if snakemake.input.get("ref"):
            samtools_opts += f" --reference {snakemake.input.ref}"

    ###################
    ### Write index ###
    ###################
    if parse_write_index:
        if "--write-index" in extra:
            sys.exit(
                "You have specified writing index (`--write-index`) in params.extra; this is automatically infered from `idx` output file."
            )
        if idx:
            samtools_opts += " --write-index"

    ###################
    ### Output file ###
    ###################
    if parse_output:
        if "-o" in extra:
            sys.exit(
                "You have specified output file (`-o`) in params.extra; this is automatically infered from the first output file."
            )
        samtools_opts += f" -o {snakemake.output[0]}"
        if idx:
            samtools_opts += f"##idx##{idx}"

    #####################
    ### Output format ###
    #####################
    if parse_output_format:
        if "-O" in extra or "--output-fmt" in extra:
            sys.exit(
                "You have specified output format (`-O/--output-fmt`) in params.extra; this is automatically infered from output file extension."
            )
        out_ext = infer_out_format(snakemake.output[0])
        samtools_opts += f" --output-fmt {out_ext}"

    return samtools_opts
