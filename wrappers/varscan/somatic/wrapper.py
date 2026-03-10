# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/09 15:42:51
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

        # Building output dirs
        makedirs(os.path.dirname(self.snakemake.output.snp))
        makedirs(os.path.dirname(self.snakemake.output.indel))

        # Output prefix
        self.prefix = os.path.splitext(self.snakemake.output.snp)[0]

        # Searching for input files
        pileup_pair = ["normal_pileup", "tumor_pileup"]

        self.in_pileup = ""
        self.mpileup = ""
        if "mpileup" in self.snakemake.input.keys():
            # Case there is a mpileup with both normal and tumor
            self.in_pileup = self.snakemake.input.mpileup
            self.mpileup = "--mpileup 1"
        elif all(pileup in self.snakemake.input.keys() for pileup in pileup_pair):
            # Case there are two separate pileup files
            self.in_pileup = " {snakemake.input.normal_pileup}" " {snakemakeinput.tumor_pileup} "
        else:
            raise KeyError("Could not find either a mpileup, or a pair of pileup files")

    def run(self):
        shell(
            "varscan somatic"  # Tool and its subcommand
            " {self.in_pileup}"  # Path to input file(s)
            " {self.prefix}"  # Path to output
            " {self.java_opts} {self.extra}"  # Extra parameters
            " {self.mpileup}"
            " --output-snp {self.snakemake.output.snp}"  # Path to snp output file
            " --output-indel {self.snakemake.output.indel}"  # Path to indel output file
        )


if __name__ == '__main__':
    Wrapper(snakemake)