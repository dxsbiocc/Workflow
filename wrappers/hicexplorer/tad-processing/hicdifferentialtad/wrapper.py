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


from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):
    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        mode = self.snakemake.params.get("mode", "all")
        assert mode in ["intra-TAD", "left-inter-TAD", "right-inter-TAD", "all"]
        if mode:
            self.extra += " --mode {}".format(mode)

        self.prefix = self.snakemake.output[0].rsplit('_')[0]

    def run(self):
        shell(
            "(hicDifferentialTAD"
            " {self.extra}"
            " --threads {self.snakemake.threads}"
            " --targetMatrix {self.snakemake.input.target}"
            " --controlMatrix {self.snakemake.input.control}"
            " --tadDomains {self.snakemake.input.domains}"
            " --outFileNamePrefix {self.prefix}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)