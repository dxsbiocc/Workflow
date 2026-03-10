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

    GENE_TABLE_DELIMITER = "\t"

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        input_args = self.snakemake.input
        if isinstance(input_args, list):
            self.input = os.path.commonpath(input_args)
        else:
            self.input = input_args

        self.output = self.snakemake.output

    def run(self):
        shell(
            "humann_join_tables"
            " --input {self.input}"
            " --output {self.output}"
            f" {self.extra}"
            " {self.log}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)
