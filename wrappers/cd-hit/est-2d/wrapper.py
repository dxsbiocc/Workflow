# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/09/04 21:54:03
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file
# ============================================================


import os
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):
    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.R1 = self.snakemake.input.get("R1")
        assert (
            isinstance(self.R1, list) and len(self.R1) == 2
        ), "R1 must be a list of two files"

        self.R1_arg = f" -i {self.R1[0]} -i2 {self.R1[1]}"

        self.R2 = self.snakemake.input.get("R2", "")
        if self.R2 and len(self.R2) == 2:
            self.R2_arg = f" -j {self.R2[0]} -j2 {self.R2[1]}"
        else:
            self.R2_arg = ""

        self.output = self.snakemake.output.get("output")
        self.paired = self.snakemake.output.get("paired")

        if self.R2 and self.paired:
            self.output_arg = f" -o {self.output} -op {self.paired}"
        else:
            self.output_arg = f" -o {self.output}"

    def run(self):
        shell("cd-hit-est-2d {self.input} {self.output} {self.extra} {self.log}")


if __name__ == "__main__":
    wrapper = Wrapper(snakemake)
