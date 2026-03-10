# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/12/09 17:11:21
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
        # Allowing for multiple FASTA files
        fasta = self.snakemake.input.get("fasta")
        assert fasta is not None, "input-> a FASTA-file is required"
        self.fasta = " ".join(fasta) if isinstance(fasta, list) else fasta


    def run(self):
        shell(
            "kallisto index"  # Tool
            " --threads {self.snakemake.threads}"
            " {self.extra}"  # Optional parameters
            " --index {self.snakemake.output.index}"  # Output file
            " {self.fasta}"  # Input FASTA files
            " {self.log}"  # Logging
        )


if __name__ == "__main__":
    Wrapper(snakemake)