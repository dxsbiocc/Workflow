# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/30 20:06:49
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from os import path
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        pass

    def run(self):
        shell(
            "(cooler digest"
            " {self.extra}"
            " -o {self.snakemake.output}"
            " {self.snakemake.input.chromsizes}"
            " {self.snakemake.input.fasta}"
            " {self.snakemake.params.enzyme}"
            ") {self.log}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)