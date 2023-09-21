# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/09/08 20:19:47
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
        self.output_format = self.snakemake.params.get("output_format", "png")
        
        background = self.snakemake.input.get("background", None)
        if background:
            self.extra = f" --backgroundModelFile {background}"
        differential = self.snakemake.input.get("differential", None)
        if differential:
            self.extra = f" --differentialTestResult {differential}"
        siginteractions = self.snakemake.input.get("siginteractions", None)
        if siginteractions:
            self.extra = f" --significantInteractions {siginteractions}"
    
    def run(self):
        shell(
            "(chicPlotViewpoint"
            " {self.extra}"
            " --threads {self.snakemake.threads}"
            " --interactionFile {self.snakemake.input.interaction}"
            " --range {self.snakemake.params.range}"
            " --outFileName {self.snakemake.output}"
            " --outputFormat {self.output_format}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)