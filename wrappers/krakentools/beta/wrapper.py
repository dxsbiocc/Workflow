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
        self.input = self.snakemake.input
        if isinstance(self.input, str):
            self.input = [self.input]
        self.output = self.snakemake.output

        # parameters
        self.type = self.snakemake.params.type
        self.level = self.snakemake.params.level

    def run(self):
        shell(
            f"beta_diversity.py"
            " {self.extra}"
            " -i {self.input}"
            " --type {self.type}"
            " --level {self.level}"
            " > {self.output}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)
