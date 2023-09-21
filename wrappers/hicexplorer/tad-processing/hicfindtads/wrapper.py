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


import os
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):
    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.correct = self.snakemake.params.get("correct", "fdr")
        assert self.correct in ["fdr", "bonferroni", "None"], "correct must be one of fdr, bonferroni, None"

        output = self.snakemake.output[0]
        if not os.path.exists(output):
            os.makedirs(output)
        prefix = self.snakemake.params.get("prefix", "")
        self.outprefix = os.path.join(output, prefix)

    def run(self):
        shell(
            "(hicFindTADs"
            " {self.extra}"
            " --numberOfProcessors {self.snakemake.threads}"
            " --correctForMultipleTesting {self.correct}"
            " --thresholdComparisons {self.snakemake.params.cutoff}"
            " --matrix {self.snakemake.input}"
            " --outPrefix {self.outprefix}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)