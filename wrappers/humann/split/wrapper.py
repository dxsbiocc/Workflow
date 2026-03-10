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
from tempfile import TemporaryDirectory
import pandas as pd
from snakemake import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.input = self.snakemake.input
        self.output_stratified = self.snakemake.output.stratified
        self.output_unstratified = self.snakemake.output.unstratified

    def run(self):
        with TemporaryDirectory() as tempdir:
            shell(
                "humann_split_stratified_table"
                " --input {self.input}"
                " --output {tempdir}"
                " {self.log}"
            )

            # mv
            shell(
                "mv {tempdir}/*_stratified* {self.output_stratified} && "
                "mv {tempdir}/*_unstratified* {self.output_unstratified}"
            )


if __name__ == "__main__":
    Wrapper(snakemake)
