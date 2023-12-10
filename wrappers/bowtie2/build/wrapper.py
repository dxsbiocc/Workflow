# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/18 19:47:36
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
        self.index = os.path.commonprefix(self.snakemake.output).rstrip(".")


    def run(self):
        shell(
            "bowtie2-build"
            " --threads {self.snakemake.threads}"
            " -f {self.snakemake.input.fasta}"
            " {self.index}"
            " {self.extra}"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)