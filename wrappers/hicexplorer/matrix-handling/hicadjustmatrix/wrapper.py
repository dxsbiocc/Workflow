# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/09/08 19:47:17
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
        regions = self.snakemake.input.get("regions", None)
        if regions:
            self.extra = f"--regions {regions}"

        chromosomes = self.snakemake.params.get("chromosomes", None)
        if chromosomes:
            self.extra += f" --chromosomes {chromosomes}"

        action = self.snakemake.params.get("action", None)
        if action:
            assert action in ['keep', 'remove', 'mask'], "action must be one of 'keep', 'remove', 'mask'"
            self.extra += f" --action {action}"
    
    def run(self):
        shell(
            "(hicAdjustMatrix"
            " {self.extra}"
            " --matrix {self.snakemake.input.matrix}"
            " --outFileName {self.snakemake.output}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)