# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/28 10:55:53
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
        cov = self.snakemake.input.get("cov", "")
        self.cov = f"-c {cov}" if cov else ""

        cutoff = self.snakemake.input.get("cutoff", "")
        self.cutoff = f"-T {cutoff}" if cutoff else ""

    def run(self):
        shell(
            "purge_dups {self.cov} {self.cutoff} {self.extra}"
            " {self.snakemake.input.paf}"
            " > {self.snakemake.output[0]}"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)
