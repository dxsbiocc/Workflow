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
from snakemake import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.input = self.snakemake.input.fastq
        # output
        self.output = self.snakemake.output
        self.outdir = os.path.dirname(self.output)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def run(self):
        with TemporaryDirectory() as tempdir:
            shell(
                "humann"
                " --input {self.input}"
                " --output {self.outdir}"
                " --threads {self.snakemake.threads}"
                " {self.extra}"
                " {self.log}"
            )


if __name__ == "__main__":
    Wrapper(snakemake)
