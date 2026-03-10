# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/25 21:37:48
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
        self.mode = self.snakemake.params.get("mode")
        assert self.mode in [
            "genome",
            "transcriptome",
            "proteins",
        ], "invalid run mode: only 'genome', 'transcriptome' or 'proteins' allowed"


        self.lineage = lineage_opt = self.snakemake.params.get("lineage", "")
        self.lineage_opt = f"--lineage {lineage_opt}" if lineage_opt else ""

        self.dataset_dir = self.snakemake.input.get("dataset_dir", "")
        if not self.dataset_dir:
            self.dataset_dir = self.snakemake.params.get("dataset_dir", "")
            if not self.dataset_dir:
                self.dataset_dir = f"dataset"

        output = os.path.abspath(self.snakemake.output[0].strip('/'))
        self.out_path = os.path.dirname(output)
        self.out = os.path.basename(output)

    def run(self):
        shell(
            "busco"
            " --cpu {self.snakemake.threads}"
            " --in {self.snakemake.input.fasta}"
            " --mode {self.mode}"
            " {self.lineage_opt}"
            " {self.extra}"
            " --download_path {self.dataset_dir}"
            " --out_path {self.out_path}"
            " --out {self.out}"
            " {self.log}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)