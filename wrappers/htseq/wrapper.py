# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2024/01/11 20:01:59
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import os
from snakemake import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        mode = self.snakemake.params.get("mode")
        if mode:
            self.extra += f" -m {mode}"
        stranded = self.snakemake.params.get("stranded")
        if stranded:
            self.extra += f" -s {stranded}"
        order = self.snakemake.params.get("order")
        if order:
            self.extra += f" -r {order}"

        self.output = os.path.splitext(self.snakemake.output.quant)[0] + ".tsv"
    
    def run(self):
        shell(
            "htseq-count --format bam"
            " -n {self.snakemake.threads}"
            " {self.extra}"
            " -c {self.output}"
            " {self.snakemake.input.bam}"
            " {self.snakemake.input.anno}"
            " {self.log}"
        )
        shell("mv {self.output} {self.snakemake.output.quant}")
    

if __name__ == '__main__':
    Wrapper(snakemake)