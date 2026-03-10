# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2024/01/11 20:01:59
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file
# ============================================================


import os
from snakemake import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.input = self.snakemake.input.kraken2_report
        self.db = self.snakemake.input.db

        self.output = self.snakemake.output.output
        self.report = self.snakemake.output.bracken_report

    def run(self):
        shell(
            "bracken"
            " {self.extra}"
            " -i {self.input}"
            " -d {self.db}"
            " -l {self.snakemake.wildcards.taxonomy}"
            " -o {self.output}"
            " -w {self.report}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)
