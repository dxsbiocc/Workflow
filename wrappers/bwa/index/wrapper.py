# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/17 17:12:46
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import os
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)
        # Check inputs/arguments.
        if len(self.snakemake.input) == 0:
            raise ValueError("A reference genome has to be provided!")
        elif len(self.snakemake.input) > 1:
            raise ValueError("Only one reference genome can be inputed!")

        # Prefix that should be used for the database
        self.prefix = self.snakemake.params.get("prefix", os.path.splitext(self.snakemake.output.idx[0])[0])

        if len(self.prefix) > 0:
            self.prefix = "-p " + self.prefix

        # Contrunction algorithm that will be used to build the database, default is bwtsw
        self.algorithm = self.snakemake.params.get("algorithm", "")

        if len(self.algorithm) != 0:
            self.algorithm = "-a " + self.algorithm

    def run(self):
        shell(
            "bwa index" 
            " {self.prefix}" 
            " {self.algorithm}" 
            " {self.snakemake.input.fasta}" 
            " {self.log}"
        )

if __name__ == '__main__':
    Wrapper(snakemake)