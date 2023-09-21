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
        plotFileName = self.snakemake.output.get("plotFileName", "p_values.txt")
        if plotFileName:
            self.extra += f" --plotFileName {plotFileName}"

        outFileName = self.snakemake.output.get("outFileName", "plot.png")
        if outFileName:
            self.extra += f" --outFileName {outFileName}"

        outFileNameData = self.snakemake.output.get("outFileNameData", "data.txt")
        if outFileNameData:
            self.extra += f" --outFileNameData {outFileNameData}"
    
    def run(self):
        shell(
            "(hicPlotSVL"
            " {self.extra}"
            " --threads {self.snakemake.threads}"
            " --matrices {self.snakemake.input}"
            " "
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)