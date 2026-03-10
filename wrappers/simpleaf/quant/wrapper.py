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
        reads1 = self.snakemake.input.fq1
        if isinstance(reads1, str):
            reads1 = [reads1]
        self.reads1 = ",".join(reads1)
        reads2 = self.snakemake.input.fq2
        if isinstance(reads2, str):
            reads2 = [reads2]
        self.reads2 = ",".join(reads2)
        self.index = self.snakemake.input.index

        resolution = self.snakemake.params.get("resolution")
        if resolution:
            self.extra += f" --resolution {resolution}"

        self.extra += f" -m {os.path.join(self.index, 't2g_3col.tsv')}"

        # options
        unfiltered = self.snakemake.input.get('unfiltered')
        if unfiltered:
            self.extra += f" --unfiltered-pl {unfiltered}"

        explicit = self.snakemake.input.get('explicit')
        if explicit:
            self.extra += f" --explicit-pl {explicit}"

        # chemistry
        self.chemistry = self.snakemake.params.get("chemistry")

    def run(self):
        path = os.path.dirname(shutil.which("simpleaf"))
        os.environ['ALEVIN_FRY_HOME'] = path
        shell(
            "simpleaf set-paths"
            " {self.log}"
        )
        if self.chemistry not in ["10xv2", "10xv3"]:
            geometry = self.snakemake.params.get("geometry")
            shell(
                "simpleaf add-chemistry"
                " --name '{self.chemistry}'"
                " --geometry '{geometry}'"
                " {self.log}"
            )
        shell(
            "(simpleaf quant "
            " --threads {self.snakemake.threads} "
            " --index {self.index}"
            " --reads1 {self.reads1}"
            " --reads2 {self.reads2}"
            " --chemistry {self.chemistry}"
            " {self.extra}"
            " --output {self.snakemake.output}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    Wrapper(snakemake)