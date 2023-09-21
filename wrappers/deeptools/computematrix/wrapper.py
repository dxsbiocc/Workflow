# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 15:17:14
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
        out_tab = self.snakemake.output.get("matrix_tab")
        out_bed = self.snakemake.output.get("matrix_bed")

        # output options
        self.optional_output = ""
        if out_tab:
            self.optional_output += " --outFileNameMatrix {out_tab} ".format(out_tab=out_tab)
        if out_bed:
            self.optional_output += " --outFileSortedRegions {out_bed} ".format(out_bed=out_bed)

        labels = self.snakemake.params.get("labels", "")
        if labels:
            labels = ' '.join([os.path.basename(label).replace(' ', '-') for label in labels])
            self.optional_output += f" --samplesLabel {labels} "

    def run(self):
        shell(
            "(computeMatrix "
            "{self.snakemake.params.command} "
            "{self.extra} "
            "-R {self.snakemake.input.bed} "
            "-S {self.snakemake.input.bigwig} "
            "-o {self.snakemake.output.matrix_gz} "
            "{self.optional_output}) {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)