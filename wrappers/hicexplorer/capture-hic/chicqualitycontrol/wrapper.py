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
        histogram = self.snakemake.output.get("histogram", None)
        if histogram:
            self.extra += " --outFileNameHistogram {}".format(histogram)
        sparsity = self.snakemake.output.get("sparsity", None)
        if sparsity:
            self.extra += " --outFileNameSparsity {}".format(sparsity)
    
    def run(self):
        shell(
            "(chicQualityControl"
            " {self.extra}"
            " --threads {self.snakemake.threads}"
            " --matrices {self.snakemake.input.matrices}"
            " --referencePoints {self.snakemake.input.reference}"
            " --sparsity {self.snakemake.params.sparsity}"
            " --outFileName {self.snakemake.output.reference}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)