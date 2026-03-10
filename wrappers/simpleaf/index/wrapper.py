# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2024/06/19 17:24:36
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import os
import shutil
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):
    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.fasta = self.snakemake.input.fasta
        self.gtf = self.snakemake.input.gtf

    def run(self):
        path = os.path.dirname(shutil.which("simpleaf"))
        os.environ['ALEVIN_FRY_HOME'] = path
        shell(
            "simpleaf set-paths"
            " {self.log}"
        )
        shell(
            "(simpleaf index "
            " --threads {self.snakemake.threads} "
            " --fasta {self.fasta}"
            " --gtf {self.gtf}"
            " {self.extra}"
            " --output {self.snakemake.output}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    Wrapper(snakemake)