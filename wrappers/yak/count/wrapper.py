# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/28 16:40:29
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import re
import os
import glob

from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.kmer = self.snakemake.params.get("kmer", 31)

        self.input_files = ""
        for file in self.snakemake.input:
            if file.endswith(".gz"):
                self.input_files += f" <(zcat {file})"
            else:
                self.input_files += f" <(cat {file})"

    def run(self):
        shell(
            "yak count"
            " {self.extra}"
            " -t{self.snakemake.threads} "
            " -k{self.kmer}"
            " -o {self.snakemake.output}"
            " {self.input_files}"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)