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
import glob
from snakemake import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        input_dir = str(self.snakemake.input)
        genefamilies = glob.glob(os.path.join(input_dir, "*2_genefamilies.tsv"))
        if len(genefamilies) == 0:
            raise ValueError("No genefamilies file found in the input directory.")
        if len(genefamilies) > 1:
            raise ValueError(
                "Multiple genefamilies files found in the input directory."
            )
        self.input = genefamilies[0]
        self.output = self.snakemake.output.ko

    def run(self):
        shell(
            "humann_regroup_table"
            " --input {self.input}"
            f" {self.extra}"
            " --output {self.output}"
            " {self.log}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)
