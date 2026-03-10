# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2024/01/16 19:59:20
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import os
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)

        self.input = str(self.snakemake.input)
        if self.input.endswith('bam'):
            self.extra += " -B"
        if self.snakemake.params.get('paired'):
            self.extra += " -P"
    
    def run(self):        
        shell(
            "(preseq lc_extrap"
            " {self.extra}"
            " -output {self.snakemake.output}"
            " {self.input}"
            ") {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)