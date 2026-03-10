# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/08 15:51:14
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from snakemake.shell import shell
from snakemake_wrapper_utils.bcftools import get_bcftools_opts
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.bcftools_opts = get_bcftools_opts(self.snakemake, parse_ref=False, parse_memory=False)

        exclude = self.snakemake.input.get("exclude", "")
        self.exclude = f"-x {exclude}" if exclude else ""
            
    def run(self):
        shell(
            "(OMP_NUM_THREADS={self.snakemake.threads} delly call"
            " -g {self.snakemake.input.ref}"
            " {self.exclude}"
            " {self.extra}"
            " {self.snakemake.input.alns} |"
            # Convert output to specified format
            " bcftools view"
            " {self.bcftools_opts}"
            ") {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)