# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/09/03 14:34:20
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
        self.log = self.snakemake.log_fmt_shell()
    
    def run(self):
        shell(
            "(wgsim {self.extra}"
            " {self.snakemake.input.ref}"
            " {self.snakemake.output.read1}"
            " {self.snakemake.output.read2}"
            ") {self.log}"
        )


if __name__ == "__main__":
    wrapper = Wrapper(snakemake)