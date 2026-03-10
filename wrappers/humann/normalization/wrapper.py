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
from tempfile import TemporaryFile
import pandas as pd
from snakemake import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.input = self.snakemake.input
        # output
        self.output = self.snakemake.output

    def run(self):
        # run humann_join_tables
        shell(
            "humann_renorm_table"
            " --input {self.input}"
            " {self.extra}"
            " --output {self.output}"
            " {self.log}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)
