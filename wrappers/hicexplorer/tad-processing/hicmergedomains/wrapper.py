# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/09/07 14:21:52
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
        protein = self.snakemake.input.get("protein", None)
        if protein:
            self.extra += " --proteinFile {}".format(protein)

        self.mergedDomains = self.snakemake.output.get("merged_domains", "mergedDomains.bed")
        self.relationList = self.snakemake.output.get("relation_list", "relationList.txt")
        self.treePlot = self.snakemake.output.get("tree_plot", "relationship_tree_")
        if not os.path.exists(self.treePlot):
            os.makedirs(self.treePlot)
        self.treePlotFormat = self.snakemake.params.get("tree_plot_format", "pdf")

    def run(self):
        shell(
            "(hicMergeDomains"
            " {self.extra}"
            " --domainFiles {self.snakemake.input}"
            " --outputMergedList {self.mergedDomains}"
            " --outputRelationList {self.relationList}"
            " --outputTreePlotPrefix {self.treePlot}/plot"
            " --outputTreePlotFormat {self.treePlotFormat}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)