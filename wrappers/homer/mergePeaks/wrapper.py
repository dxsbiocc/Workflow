# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/21 15:03:24
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        if "-prefix" in self.extra:
            raise ValueError(
                "The use of the -prefix parameter is not yet supported in this wrapper"
            )
            
    def run(self):
        shell(
            "(mergePeaks" 
            " {self.snakemake.input}" 
            " {self.extra}" 
            " > {self.snakemake.output})" 
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)