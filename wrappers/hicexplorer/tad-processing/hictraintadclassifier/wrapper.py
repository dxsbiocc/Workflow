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
        self.mode = self.snakemake.params.get("mode")
        assert self.mode in [
            "train_new", "train_existing", 
            "train_test", "predict_test"
        ], "mode must be in ['train_new', 'train_existing', 'train_test', 'predict_test']"
        self.normalization_method = self.snakemake.params.get("normalization_method")
        assert self.normalization_method in ["obs_exp", "range"], "normalization_method must be in ['obs_exp', 'range']"
        self.resolution = self.snakemake.params.get("resolution")

        saved_classifier = self.snakemake.output.get("classifier")
        if saved_classifier:
            self.extra = f" --saved_classifier {saved_classifier}"
        protein_file = self.snakemake.input.get("protein_file")
        if protein_file:
            self.extra = f" --protein_file {protein_file}"

    def run(self):
        shell(
            "(hicTrainTADClassifier"
            " {self.extra}"
            " --threads {self.snakemake.threads}"
            " --mode {self.mode}"
            " --normalization_method {self.normalization_method}"
            " --resolution {self.resolution}"
            " --matrices {self.snakemake.input.matrices}"
            " --domain_file {self.snakemake.input.domain}"
            " --out_file {self.snakemake.output.out}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)