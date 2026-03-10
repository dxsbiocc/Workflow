# -*- encoding: utf-8 -*-
# ============================================================
# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/09/04 15:09:36
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
        rangeInBase = self.snakemake.params.get("rangeInBase", "")
        rangeInBins = self.snakemake.params.get("rangeInBins", "")
        assert not (rangeInBase and rangeInBins), "Only one of --range and --rangeInBins can be specified."

        if rangeInBase:
            self.extra += f" --range {rangeInBase}"

        if rangeInBins:
            self.extra += f" --rangeInBins {rangeInBins}"
    
    def run(self):
        shell(
            "(hicAverageRegions"
            " --matrix {self.snakemake.input.mat}"
            " --regions {self.snakemake.input.bed}"
            " {self.extra}"
            " --outFileName {self.snakemake.output}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)