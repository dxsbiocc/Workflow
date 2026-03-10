# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2024/01/11 20:01:59
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file
# ============================================================


import os
from tempfile import TemporaryDirectory
from snakemake import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.input = self.snakemake.input.fastq
        # output
        self.outdir = self.snakemake.output
        if not os.path.exists(str(self.outdir)):
            os.makedirs(str(self.outdir))

        self.nucleotide = self.snakemake.params.get("nucleotide", "")
        self.protein = self.snakemake.params.get("protein", "")
        self.utility_mapping = self.snakemake.params.get("utility_mapping", "")

        self.extra = self.snakemake.params.get("extra", "")
        metaphlan = self.snakemake.params.get("metaphlan", "")
        if metaphlan:
            self.extra += f" --metaphlan-options '{metaphlan}'"

    def run(self):
        # update the nucleotide, protein, utility_mapping
        shell("humann_config --update database_folders nucleotide {self.nucleotide}")
        shell("humann_config --update database_folders protein {self.protein}")
        shell(
            "humann_config --update database_folders utility_mapping {self.utility_mapping}"
        )
        # run humann
        shell(
            "humann"
            " --input {self.input}"
            " --output {self.outdir}"
            " --threads {self.snakemake.threads}"
            " {self.extra}"
            " --o-log {self.snakemake.log}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)
