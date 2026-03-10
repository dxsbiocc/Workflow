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
        self.output = self.snakemake.output
        # parameters
        self.subcommand = self.snakemake.params.get("subcommand")
        if self.subcommand not in ["kreport2mpa", "kreport2krona"]:
            raise ValueError("subcommand must be kreport2mpa or kreport2krona")

    def run(self):
        shell(
            f"{self.subcommand}.py"
            " {self.extra}"
            " -r {self.input}"
            " -o {self.output}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)
