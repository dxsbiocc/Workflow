# -*- encoding: utf-8 -*-
# ============================================================
# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/09/04 15:09:36
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
        chrPrefixLoops = self.snakemake.params.get('', '')
        if chrPrefixLoops:
            self.extra += f" --chrPrefixLoops {chrPrefixLoops}"
    
    def run(self):
        shell(
            "(hicHyperoptDetectLoops"
            " {self.extra}"
            " --matrix {self.snakemake.input.mat}"
            " --proteinFile {self.snakemake.input.protein}"
            " --maximumNumberOfLoops {self.snakemake.params.ml}"
            " --outputFileName {self.snakemake.output}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)