# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/10/20 17:24:26
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
        out_ext = self.snakemake.output[0].split('.')[-1]
        if out_ext in ['bw', 'bigwig', 'bigWig']:
            self.out_format = 'bigwig'
        elif out_ext in ['bedgraph', 'bedGraph']:
            self.out_format = 'bedgraph'
        else:
             raise TypeError('Output extension should be one of ["bw", "bigwig", "bigWig", "bedgraph", "bedGraph"]')

    def run(self):
        shell(
            "( bamCoverage "
            "--bam {self.snakemake.input} "
            "--outFileName {self.snakemake.output} "
            "--outFileFormat {self.out_format} "
            "--numberOfProcessors {self.snakemake.threads} "
            "{self.extra} "
            ") "
            "{log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)