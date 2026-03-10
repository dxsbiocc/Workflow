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
        bw = self.snakemake.input.get('bigwig')
        if bw:
            self.extra += f" --bigwig {bw}"
        loops = self.snakemake.input.get('loops')
        if loops:
            self.extra += f" --loops {loops}"
        tads = self.snakemake.input.get('tads')
        if tads:
            self.extra += f" --tads {tads}"

        title = self.snakemake.params.get('title')
        if title:
            self.extra += f" --title {title}"
        
    
    def run(self):
        shell(
            "(hicPlotMatrix"
            " {self.extra}"
            " --matrix {self.snakemake.input.matrix}"
            " --outFileName {self.snakemake.output}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)