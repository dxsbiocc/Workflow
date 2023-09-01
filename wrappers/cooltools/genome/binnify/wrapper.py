# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/31 17:28:07
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
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)
        ## Extract arguments
        self.binsize = self.snakemake.params.get("binsize", self.snakemake.wildcards.get("binsize", 0))
        if not self.binsize:
            raise ValueError("Please specify binsize either as a wildcard or as a parameter")

    def run(self):
        shell(
            "(cooltools genome binnify"
            " {self.snakemake.input.chromsizes}"
            " {self.binsize} "
            " {self.extra} "
            " > {self.snakemake.output})"
            " {self.log}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)