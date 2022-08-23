import sys


def get_java_opts(snakemake):
    """Obtain java_opts from params, and handle resource definitions in resources."""

    java_opts = snakemake.params.get("java_opts", "")
    extra = snakemake.params.get("extra", "")

    # Getting memory in megabytes, if java opts is not filled with -Xmx parameter
    # By doing so, backward compatibility is preserved
    if "mem_mb" in snakemake.resources.keys():
        if "-Xmx" in java_opts:
            sys.exit(
                "You have specified resources.mem_mb and provided `-Xmx` in params.java_opts. For Java memory specifications, please only use resources.mem_mb."
            )
        if "-Xmx" in extra:
            sys.exit(
                "You have specified resources.mem_mb and provided `-Xmx` in params.extra. For Java memory specifications, please only use resources.mem_mb."
            )
        java_opts += " -Xmx{}M".format(snakemake.resources["mem_mb"])

    # Getting memory in gigabytes, for user convenience. Please prefer the use
    # of mem_mb over mem_gb as advised in documentation.
    elif "mem_gb" in snakemake.resources.keys():
        if "-Xmx" in java_opts:
            sys.exit(
                "You have specified resources.mem_gb and provided `-Xmx` in params.java_opts. For Java memory specifications, please only use resources.mem_mb."
            )
        if "-Xmx" in extra:
            sys.exit(
                "You have specified resources.mem_gb and provided `-Xmx` in params.extra. For Java memory specifications, please only use resources.mem_mb."
            )
        java_opts += " -Xmx{}G".format(snakemake.resources["mem_gb"])

    # Getting java temp directory from output files list, if -Djava.io.tmpdir
    # is not provided in java parameters. By doing so, backward compatibility is
    # not broken.
    if "java_temp" in snakemake.output.keys():
        if "-Djava.io.tmpdir" in java_opts:
            sys.exit(
                "You have specified output.java_temp and provided `-Djava.io.tmpdir` in params.java_opts. Please choose the one you intended and remove the other specification."
            )
        if "-Djava.io.tmpdir" in extra:
            sys.exit(
                "You have specified output.java_temp and provided `-Djava.io.tmpdir` in params.extra. Please choose the one you intended and remove the other specification."
            )
        java_opts += " -Djava.io.tmpdir={}".format(snakemake.output["java_temp"])

    return java_opts
