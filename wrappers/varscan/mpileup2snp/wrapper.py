# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/09 15:56:01
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import os
from snakemake.shell import shell
from snakemake.utils import makedirs
from snakemake_wrapper_utils.java import get_java_opts
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.java_opts = get_java_opts(self.snakemake)
        # In case input files are gzipped mpileup files,
        # they are being unzipped and piped
        # In that case, it is recommended to use at least 2 threads:
        # - One for unzipping with zcat
        # - One for running varscan
        self.pileup = (
            " cat {} ".format(self.snakemake.input[0])
            if not self.snakemake.input[0].endswith("gz")
            else " zcat {} ".format(self.snakemake.input[0])
        )

        # Building output directories
        makedirs(os.path.dirname(self.snakemake.output[0]))

    def run(self):
        shell(
            "varscan mpileup2snp "  # Tool and its subprocess
            "<( {self.pileup} ) "
            "{self.java_opts} {self.extra} "  # Extra parameters
            "> {self.snakemake.output[0]} "  # Path to vcf file
            "{self.log}"  # Logging behaviour
        )


if __name__ == '__main__':
    Wrapper(snakemake)