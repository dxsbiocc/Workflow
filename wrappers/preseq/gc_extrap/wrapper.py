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
        max_width = self.input.endswith('max_width')
        if max_width:
            self.extra += " -max_width {max_width}"
        bin_size = self.snakemake.params.get('bin_size')
        if bin_size:
            self.extra += " -bin_size {bin_size}"
    
    def run(self):        
        shell(
            "(preseq gc_extrap"
            " {self.extra}"
            " -output {self.snakemake.output}"
            " {self.input}"
            ") {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)