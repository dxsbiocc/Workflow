# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/10/20 15:56:35
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
        pass

    def run(self):
        # index
        if not os.path.exists(self.snakemake.input[0] + '.bai'):
            shell("samtools index {self.snakemake.input.bam} ")

        shell(
            "(alignmentSieve "
            "--numberOfProcessors {self.snakemake.threads} "
            "--bam {self.snakemake.input[0]} "
            "-o {self.snakemake.output} "
            "{self.extra}) {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)