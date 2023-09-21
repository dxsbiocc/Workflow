# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 11:27:09
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import os
from snakemake.shell import shell
from snakemake_wrapper_utils.samtools import get_samtools_opts
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)

        # gc file
        gc = self.snakemake.input.get("gc", "")
        self.gc = f"-r {gc}" if gc else ""

    def run(self):
        shell(
            "plot-bamstats {self.extra} "
            "{self.gc} "
            "-p {self.snakemake.output}/ "
            "{self.snakemake.input.stats} "
            "{self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)