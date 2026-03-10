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
        self.normalization_method = self.snakemake.params.get("normalization_method")
        assert self.normalization_method in ["obs_exp", "range"], "normalization_method must be in ['obs_exp', 'range']"

        saved_classifier = self.snakemake.input.get("classifier")
        if saved_classifier:
            self.extra = f" --saved_classifier {saved_classifier}"

    def run(self):
        shell(
            "(hicTADClassifier"
            " {self.extra}"
            " --threads {self.snakemake.threads}"
            " --normalization_method {self.normalization_method}"
            " --matrices {self.snakemake.input.matrices}"
            " --out_file {self.snakemake.output.out}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)