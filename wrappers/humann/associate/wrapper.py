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
        self.input = self.snakemake.input.abundance
        self.group_line = self.snakemake.params.group
        # output
        self.output = self.snakemake.output.associate
        self.pcl = self.snakemake.output.pcl

    def run(self):
        # generate pcl file
        with open(str(self.input), "r") as f:
            lines = f.readlines()
        lines.insert(1, f"Group\t{self.group_line}\n")
        with open(self.pcl, "w") as f:
            f.writelines(lines)
        # run humann_associate
        shell(
            "humann_associate"
            " --input {self.pcl}"
            " --focal-metadatum Group --focal-type categorical"
            " --last-metadatum Group --fdr 0.2"
            " --output {self.output}"
            " {self.log}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)
