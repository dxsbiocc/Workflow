# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/10 10:07:50
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from os.path import commonprefix
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake):
        super().__init__(snakemake)

    def parser(self):
        pass

    def run(self):
        # Detemining common prefix in output files
        # to fill the requested parameter '-o'
        prefix = commonprefix(self.snakemake.output)

        shell(
            "msisensor msi"  # Tool and its sub-command
            " -d {self.snakemake.input.microsat}"  # Path to homopolymer/microsat file
            " -n {self.snakemake.input.normal}"  # Path to normal bam
            " -t {self.snakemake.input.tumor}"  # Path to tumor bam
            " -o {prefix}"  # Path to output distribution file
            " -b {self.snakemake.threads}"  # Maximum number of threads used
            " {self.extra}"  # Optional extra parameters
            " {self.log}"  # Logging behavior
        )


if __name__ == '__main__':
    Wrapper(snakemake)