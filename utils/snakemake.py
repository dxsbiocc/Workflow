def get_mem(snakemake, out_unit="MiB"):
    """
    Obtain requested memory (from resources) and return in given units.
    If no memory resources found, return 0.
    """

    # Store memory in MiB
    mem_mb = snakemake.resources.get("mem_gb", 0) * 1024
    if not mem_mb:
        mem_mb = snakemake.resources.get("mem_mb", 0)

    if out_unit == "KiB":
        return mem_mb * 1024
    elif out_unit == "MiB":
        return mem_mb
    elif out_unit == "GiB":
        return mem_mb / 1024
    else:
        raise ValueError("invalid output unit. Only KiB, MiB and GiB supported.")
