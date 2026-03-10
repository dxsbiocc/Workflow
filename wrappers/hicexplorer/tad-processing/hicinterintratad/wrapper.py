# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/09/07 14:21:52
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):
    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.outFileNameRatioPlot = self.snakemake.output.get("ratioplot", "ratio.png")

    def run(self):
        shell(
            "(hicInterIntraTAD"
            " {self.extra}"
            " --threads {self.snakemake.threads}"
            " --matrix {self.snakemake.input.matrix}"
            " --tadDomains {self.snakemake.input.domains}"
            " --outFileName {self.snakemake.output.out}"
            " --outFileNameRatioPlot {self.outFileNameRatioPlot}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)