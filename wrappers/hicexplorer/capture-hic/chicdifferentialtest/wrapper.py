# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/09/08 20:19:47
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
        self.test = self.snakemake.params.get("test", "fisher")
        assert self.test in ["fisher", "chi2"], "test must be one of fisher, chi2"
    
    def run(self):
        shell(
            "(chicDifferentialTest"
            " {self.extra}"
            " --threads {self.snakemake.threads}"
            " --aggregatedFile {self.snakemake.input}"
            " --alpha {self.snakemake.params.alpha}"
            " --statisticTest {self.test}"
            " --outFileName {self.snakemake.output}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)