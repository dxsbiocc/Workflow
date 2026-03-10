# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/31 17:31:41
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
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)

    def run(self):
        shell(
            "(cooltools genome gc"
            " {self.extra}"
            " {self.snakemake.input.bins}"
            " {self.snakemake.input.fasta}"
            " > {self.snakemake.output})"
            " {self.log} "
        )


if __name__ == "__main__":
    Wrapper(snakemake)