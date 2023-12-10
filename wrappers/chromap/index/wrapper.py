# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/09/21 14:50:15
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
        kmer = self.snakemake.params.get("kmer", "")
        self.extra += f"-k {kmer}" if kmer else ""

        mwsize = self.snakemake.params.get("miniWinSize", "")
        self.extra += f"-w {mwsize}" if mwsize else ""

    def run(self):
        shell(
            "chromap -i"
            " {self.extra}"
            " -r {self.snakemake.input}"
            " -o {self.snakemake.output}"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)