# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 15:09:00
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
            self.snakemake, parse_write_index=False, parse_output_format=False
        )

        # check for optional bed file
        bed = self.snakemake.input.get("bed", None)
        self.bed = "-b " + bed if bed else ""

    def run(self):
        shell("samtools depth"
        " {self.samtools_opts} {self.extra}"
        " {self.bed}"
        " {self.snakemake.input.bams}"
        " {self.log}")


if __name__ == '__main__':
    Wrapper(snakemake)