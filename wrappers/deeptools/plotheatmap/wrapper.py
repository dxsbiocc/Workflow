# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 15:46:04
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
        out_region = self.snakemake.output.get("regions")
        out_matrix = self.snakemake.output.get("heatmap_matrix")

        self.optional_output = ""
        if out_region:
            self.optional_output += f" --outFileSortedRegions {out_region} "
        if out_matrix:
            self.optional_output += f" --outFileNameMatrix {out_matrix} "
        labels = self.snakemake.params.get("labels", "")
        if labels:
            labels = ' '.join([os.path.basename(label).replace(' ', '-') for label in labels])
            self.optional_output += f" --samplesLabel {labels} "

    def run(self):
        shell(
            "(plotHeatmap "
            "-m {self.snakemake.input[0]} "
            "-o {self.snakemake.output.heatmap_img} "
            "{self.optional_output} "
            "{self.extra}) {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)