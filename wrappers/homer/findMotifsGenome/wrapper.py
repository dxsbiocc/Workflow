# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/22 10:52:09
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
            "(findMotifsGenome.pl"
            " {self.snakemake.input.peak}"
            " {self.snakemake.input.genome}"
            " {self.snakemake.output}"
            " {self.extra}"
            " -p {self.snakemake.threads}"
            ") {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)