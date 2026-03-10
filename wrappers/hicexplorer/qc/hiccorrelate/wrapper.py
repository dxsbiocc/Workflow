# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/31 20:41:02
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
        self.method = self.snakemake.params.get("method", "pearson")
        assert self.method in ["pearson", "spearman"], "method must be 'pearson' or 'spearman'"
        
        assert len(self.snakemake.input) == len(self.snakemake.params.labels), \
            "labels must be the same length as input matrices"
    
    def run(self):
        shell(
            "(hicCorrelate"
            " {self.extra}"
            " --matrices {self.snakemake.input}"
            " --method={self.method}"
            " --labels {self.snakemake.params.labels}"
            " --colorMap {self.snakemake.params.colorMap}"
            " --outFileNameHeatmap {self.snakemake.output.heatmap}"
            " --outFileNameScatter {self.snakemake.output.scatter}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)