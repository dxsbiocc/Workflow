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
        return super().parser()
    
    def run(self):
        shell(
            "(chicSignificantInteractions"
            " {self.extra}"
            " --interactionFile {self.snakemake.input.interaction}"
            " --backgroundModelFile {self.snakemake.input.background}"
            " --range {self.snakemake.params.range}"
            " --pValue {self.snakemake.params.pvalue}"
            " --outFileNameSignificant {self.snakemake.output.significant}"
            " --outFileNameTarget {self.snakemake.output.target}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)