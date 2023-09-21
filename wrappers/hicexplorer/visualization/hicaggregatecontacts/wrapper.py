# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/09/08 14:47:04
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
        self.mode = self.snakemake.params.get('mode')
        assert self.mode in ['inter-chr', 'intra-chr', 'all'], "mode must be one of ['inter-chr', 'intra-chr', 'all']"
        
        heapmap = self.snakemake.output.get("heatmap")
        if heapmap:
            self.extra += " --diagnosticHeatmapFile {heapmap}".format(heapmap=heapmap)

        premat = self.snakemake.params.get("premat")
        if premat:
            self.extra += " --outFilePrefixMatrix {premat}".format(premat=premat)
        contact = self.snakemake.params.get("contact")
        if contact:
            self.extra += " --outFileContactPairs {contact}".format(contact=contact)
    
    def run(self):
        shell(
            "(hicAggregateContacts"
            " {self.extra}"
            " --matrix {self.snakemake.input.matrix}"
            " --BED {self.snakemake.input.bed}"
            " --mode {self.mode}"
            " --outFileName {self.snakemake.output.agg}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)