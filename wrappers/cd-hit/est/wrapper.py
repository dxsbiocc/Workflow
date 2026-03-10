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
        self.input = self.snakemake.input
        if isinstance(self.input, list) and len(self.input) == 2:
            self.input = f" -i {self.input[0]} -j {self.input[1]}"
        else:
            self.input = f" -i {self.input}"

        self.output = self.snakemake.output
        if isinstance(self.output, list) and len(self.output) == 2:
            self.output = f" -o {self.output[0]} -op {self.output[1]}"
        else:
            self.output = f" -o {self.output}"

    def run(self):
        shell("cd-hit-est {self.input} {self.output} {self.extra} {self.log}")


if __name__ == "__main__":
    wrapper = Wrapper(snakemake)
