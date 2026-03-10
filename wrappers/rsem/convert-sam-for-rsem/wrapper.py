# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/04/21 14:12:19
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
        self.input_bam = self.snakemake.input.get("bam")
        self.output = os.path.splitext(self.snakemake.output.get("bam"))[0]

    def run(self):
        shell(
            "convert-sam-for-rsem"
            " --num-threads {self.snakemake.threads}"
            " {self.input_bam}"
            " {self.output}"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)