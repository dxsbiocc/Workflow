# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/08 21:32:00
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import os
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        # Allowing for multiple FASTA files
        fasta = self.snakemake.input.get("fasta")
        assert fasta is not None, "input-> a FASTA-file or a sequence is required"
        self.input_seq = ""
        if not "." in fasta:
            self.input_seq += "-c "
        self.input_seq += ",".join(fasta) if isinstance(fasta, list) else fasta

        self.output = self.snakemake.output[0]
        if not os.path.exists(self.output):
            os.makedirs(self.output)
        self.prefix = os.path.join(self.output, self.snakemake.params.get("prefix", "hisat2_index"))
        

    def run(self):
        shell(
            "hisat2-build {self.extra}"
            " -p {self.snakemake.threads}"
            " {self.input_seq}"
            " {self.prefix}"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)