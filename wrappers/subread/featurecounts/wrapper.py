# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/21 15:07:19
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import os
import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        # optional input files and directories
        self.strand = self.snakemake.params.get("strand", 0)

        self.fasta = self.snakemake.input.get("fasta", "")
        if self.fasta:
            self.fasta = f"-G {self.fasta}"

        self.chr_names = self.snakemake.input.get("chr_names", "")
        if self.chr_names:
            self.chr_names = f"-A {self.chr_names}"

        self.r_path = self.snakemake.params.get("r_path", "")
        if self.r_path:
            self.r_path = f"--Rpath {self.r_path}"

        paired = self.snakemake.params.get("paired")
        if paired:
            self.extra += " -p --countReadPairs"

        self.out = self.snakemake.output.get("quant")

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "featureCounts"
                " -T {self.snakemake.threads}"
                " -s {self.strand}"
                " -a {self.snakemake.input.annotation}"
                " {self.fasta}"
                " {self.chr_names}"
                " {self.r_path}"
                " {self.extra}"
                " --tmpDir {tmpdir}"
                " -o {self.out}"
                " {self.snakemake.input.samples}"
                " {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)