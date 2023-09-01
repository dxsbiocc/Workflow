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
        self.command = self.snakemake.params.get("command", "")
        assert self.command in [
            "diagnostic_plot", 
            "correct"
        ], "`command` must be 'diagnostic_plot' or 'correct'"

        if self.command == "correct":
            correctionMethod = self.snakemake.params.get("correctionMethod", "KR")
            self.extra += f" --correctionMethod {correctionMethod}"
            if correctionMethod == "ICE":
                ICE = self.snakemake.params.get("ICE", {})
                for k, v in ICE:
                    self.extra += f"{k} {v}"

        chromosomes = self.snakemake.params.get("chromosomes", [])
        if chromosomes:
            self.extra += " --chromosomes " + " ".join(chromosomes)

        perchr = self.snakemake.params.get("perchr", False)
        if perchr:
            self.extra += " --perchr"
    
    def run(self):
        shell(
            "(hicCorrectMatrix {self.command}"
            " {self.extra}"
            " --matrix {self.snakemake.input}"
            " -o {self.snakemake.output}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)