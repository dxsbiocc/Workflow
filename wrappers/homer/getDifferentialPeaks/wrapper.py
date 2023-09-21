# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/21 14:59:45
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
        return super().parser()
    
    def run(self):
        shell(
            "(getDifferentialPeaks"
            " {self.snakemake.input.peaks}"
            " {self.snakemake.input.first}"
            " {self.snakemake.input.second}"
            " {self.extra}"
            " > {self.snakemake.output})"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)