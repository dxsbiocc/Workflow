# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 14:57:22
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


from snakemake.shell import shell
from snakemake_wrapper_utils.samtools import get_samtools_opts
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.samtools_opts = get_samtools_opts(
            self.snakemake, parse_threads=False, parse_write_index=False, parse_output_format=False
        )

    def run(self):
        shell("samtools faidx"
        " {self.samtools_opts} {self.extra}"
        " {self.snakemake.input[0]} {self.log}")


if __name__ == '__main__':
    Wrapper(snakemake)