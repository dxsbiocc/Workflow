# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/25 16:22:12
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import os
import sys
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.kmer = self.snakemake.params.get("kmer", 31)        

    def run(self):
        shell(
            "(meryl k={self.kmer} count"
            " {self.snakemake.input}"
            " output {self.snakemake.output}"
            " {self.extra}"
            " ) {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)