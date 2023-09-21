# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/31 20:41:02
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
        assert len(self.snakemake.input) == len(self.snakemake.params.labels), \
            "labels must be the same length as input matrices"
    
    def run(self):
        shell(
            "(hicPlotDistVsCounts"
            " {self.extra}"
            " --matrices {self.snakemake.input}"
            " --labels {self.snakemake.params.labels}"
            " --plotFile {self.snakemake.output}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)