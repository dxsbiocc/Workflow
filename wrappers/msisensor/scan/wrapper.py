# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/10 10:03:48
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

    def __init__(self, snakemake):
        super().__init__(snakemake)

    def parser(self):
        pass

    def run(self):
        shell(
            "msisensor scan "  # Tool and its sub-command
            "-d {self.snakemake.input} "  # Path to fasta file
            "-o {self.snakemake.output} "  # Path to output file
            "{self.extra} "  # Optional extra parameters
            "{self.log}"  # Logging behavior
        )


if __name__ == '__main__':
    Wrapper(snakemake)