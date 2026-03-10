# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 09:29:47
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file
# ============================================================


import os
import sys
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.input = self.snakemake.input
        self.output = self.snakemake.output

    def run(self):
        shell("lefse_run.py {self.extra} {self.input} {self.output}")


if __name__ == "__main__":
    Wrapper(snakemake)
