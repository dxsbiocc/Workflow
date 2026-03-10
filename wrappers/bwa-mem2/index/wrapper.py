# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/17 17:12:46
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import os
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        # Check inputs/arguments.
        self.prefix = os.path.commonprefix(self.snakemake.output.get('index'))
        dirname = os.path.dirname(self.prefix)
        if not os.path.exists(dirname):
            os.makedirs(dirname)

    def run(self):
        shell(
            "bwa-mem2 index" 
            " -p {self.prefix}" 
            " {self.snakemake.input.fasta}" 
            " {self.log}"
        )

if __name__ == '__main__':
    Wrapper(snakemake)