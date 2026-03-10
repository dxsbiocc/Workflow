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
        self.input = self.snakemake.input.fastq
        self.db = self.snakemake.input.db
        self.output = self.snakemake.output.output
        self.report = self.snakemake.output.report
        self.extra = self.snakemake.params.extra

    def run(self):
        shell(
            "kraken2"
            " {self.extra} "
            " --paired {self.input}"
            " --threads {self.snakemake.threads}"
            " --db {self.db}"
            " --report {self.report}"
            " --output {self.output}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)
