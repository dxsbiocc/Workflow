# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 19:57:37
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import os
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)

        self.genome = ""
        self.input_file = ""

        bam = self.snakemake.input.get("bam", None)
        bed = self.snakemake.input.get("bed", None)
        if bool(bam) ^ bool(bed): 
            if bam:
                self.input_file = f"-ibam {bam}"
            else:
                genome = self.snakemake.input.get("ref", None)
                assert genome, "must have genome file!"
                self.input_file = f"-i {bed}"
                self.genome = f"-g {genome}"
        else:
            raise ValueError("bam or bed must be choose one!")


    def run(self):
        shell(
            "(genomeCoverageBed"
            " {self.input_file}"
            " {self.genome}"
            " {self.extra}"
            " > {self.snakemake.output[0]}) {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)