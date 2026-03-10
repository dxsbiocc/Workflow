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
        pass
    
    def run(self):
        shell(
            "(hicDetectLoops"
            " --matrix {self.snakemake.input}"
            " --outFileName {self.snakemake.output}"
            " {self.extra}"
            " --threads {self.snakemake.threads}"
            " --peakWidth {self.snakemake.params.peakWidth}"
            " --windowSize {self.snakemake.params.windowSize}"
            " --expected {self.snakemake.params.expected}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)