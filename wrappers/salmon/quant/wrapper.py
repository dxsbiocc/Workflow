# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/12/09 16:39:35
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from os.path import dirname
from snakemake.shell import shell

from snakemake_wrapper_utils.base import WrapperBase


class MixedPairedUnpairedInput(Exception):
    def __init__(self):
        super().__init__(
            "Salmon cannot quantify mixed paired/unpaired input files. "
            "Please input either `r1`, `r2` (paired) or `r` (unpaired)"
        )


class MissingMateError(Exception):
    def __init__(self):
        super().__init__(
            "Salmon requires an equal number of paired reads in `r1` and `r2`,"
            " or a list of unpaired reads `r`"
        )


def uncompress_bz2(snake_io, salmon_threads):
    """
    Provide bzip2 on-the-fly decompression

    For each of these b-unzipping, a thread will be used. Therefore, the maximum number of threads given to Salmon
    shall be reduced by one in order not to be killed on a cluster.
    """

    # Asking forgiveness instead of permission
    try:
        # If no error are raised, then we have a string.
        if snake_io.endswith("bz2"):
            return [f"<( bzip2 --decompress --stdout {snake_io} )"], salmon_threads - 1
        return [snake_io], salmon_threads
    except AttributeError:
        # As an error has been raise, we have a list of fastq files.
        fq_files = []
        for fastq in snake_io:
            if fastq.endswith("bz2"):
                fq_files.append(f"<( bzip2 --decompress --stdout {fastq} )")
                salmon_threads -= 1
            else:
                fq_files.append(fastq)
        return fq_files, salmon_threads


class Wrapper(WrapperBase):
    def __init__(self, snakemake):
        super().__init__(snakemake)

    def parser(self):
        self.libtype = self.snakemake.params.get("libtype", "A")
        self.max_threads = self.snakemake.threads

        if "--validateMappings" in self.extra:
            raise DeprecationWarning("`--validateMappings` is deprecated and has no effect")

        r1 = self.snakemake.input.get("r1")
        r2 = self.snakemake.input.get("r2")
        r = self.snakemake.input.get("r")

        if all(mate is not None for mate in [r1, r2]):
            r1, self.max_threads = uncompress_bz2(r1, self.max_threads)
            r2, self.max_threads = uncompress_bz2(r2, self.max_threads)

            if len(r1) != len(r2):
                raise MissingMateError()
            if r is not None:
                raise MixedPairedUnpairedInput()

            r1_cmd = " --mates1 {}".format(" ".join(r1))
            r2_cmd = " --mates2 {}".format(" ".join(r2))
            self.read_cmd = " ".join([r1_cmd, r2_cmd])

        elif r is not None:
            if any(mate is not None for mate in [r1, r2]):
                raise MixedPairedUnpairedInput()

            r, self.max_threads = uncompress_bz2(r, self.max_threads)
            self.read_cmd = " --unmatedReads {}".format(" ".join(r))

        else:
            raise MissingMateError()

        self.gene_map = self.snakemake.input.get("gtf", "")
        if self.gene_map:
            self.gene_map = f"--geneMap {self.gene_map}"

        self.bam = self.snakemake.output.get("bam", "")
        if self.bam:
            self.bam = f"--writeMappings {self.bam}"

        self.outdir = self.snakemake.output.get("quant")
        self.index = self.snakemake.input["index"]


        if isinstance(self.index, list):
            self.index = dirname(self.index[0])

        if self.max_threads < 1:
            raise ValueError(
                "On-the-fly b-unzipping have raised the required number of threads. "
                f"Please request at least {1 - self.max_threads} more threads."
            )
        
    def run(self):
        shell(
            "salmon quant --index {self.index} "
            " --libType {self.libtype} {self.read_cmd} --output {self.outdir} {self.gene_map} "
            " --threads {self.max_threads} {self.extra} {self.bam} {self.log}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)