# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 11:43:23
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
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)

        if not self.snakemake.output[0].endswith(".gz"):
            raise Exception(
                'output file will be compressed and therefore filename should end with ".gz"'
            )

    def run(self):
        shell(
            "(samtools mpileup {self.extra} -f {self.snakemake.input.reference} {self.snakemake.input.bam} \
                | pigz > {self.snakemake.output}) {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)