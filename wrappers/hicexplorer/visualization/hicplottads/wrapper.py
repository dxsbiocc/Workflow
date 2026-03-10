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
        bed = self.snakemake.input.get("bed", None)
        if bed:
            self.extra += f" --BED {bed}"

        title = self.snakemake.params.get("title", None)
        if title:
            self.extra += f" --title {title}"
    
    def run(self):
        shell(
            "(hicPlotTADs"
            " {self.extra}"
            " --tracks {self.snakemake.input.tracks}"
            " --outFileName {self.snakemake.output}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)