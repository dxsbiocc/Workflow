# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 21:34:12
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):

        self.cat = None
        if len(self.snakemake.input) > 1:
            if all(f.endswith(".gz") for f in self.snakemake.input):
                self.cat = "zcat"
            elif all(not f.endswith(".gz") for f in self.snakemake.input):
                self.cat = "cat"
            else:
                raise ValueError("Input files must be all compressed or uncompressed.")

    def run(self):
        if self.cat:
            shell(
                "({self.cat} {self.snakemake.input} | "
                "sort -k1,1 -k2,2n | "
                "bedtools merge {self.extra} "
                "-i stdin > {self.snakemake.output}) "
                " {self.og}"
            )
        else:
            shell(
                "( bedtools merge"
                " {self.extra}"
                " -i {self.snakemake.input}"
                " > {self.snakemake.output})"
                " {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)