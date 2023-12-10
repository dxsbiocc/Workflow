# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/12/09 17:14:45
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from snakemake.shell import shell

from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        # Allowing for multiple FASTQ files
        fastq = self.snakemake.input.get("fastq")
        assert fastq is not None, "input-> a FASTQ-file is required"
        self.fastq = " ".join(fastq) if isinstance(fastq, list) else fastq

    def run(self):
        shell(
            "kallisto quant "  # Tool
            " --threads {self.snakemake.threads}"  # Number of threads
            " --index {self.snakemake.input.index}"  # Input file
            " {self.extra}"  # Optional parameters
            " --output-dir {self.snakemake.output}"  # Output directory
            " {self.fastq}"  # Input FASTQ files
            " {self.log}"  # Logging
        )


if __name__ == "__main__":
    Wrapper(snakemake)