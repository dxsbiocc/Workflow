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


import tempfile
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

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dataset_dir = self.snakemake.input.get("dataset_dir", "")
            if not dataset_dir:
                dataset_dir = f"{tmpdir}/dataset"

            shell(
                "busco"
                " --cpu {self.snakemake.threads}"
                " --in {self.snakemake.input}"
                " --mode {self.mode}"
                " {self.lineage_opt}"
                " {self.extra}"
                " --download_path {dataset_dir}"
                " --out_path {tmpdir}"
                " --out output"
                " {self.log}"
            )

            if self.snakemake.output.get("short_txt"):
                assert self.lineage, "parameter 'lineage' is required to output 'short_tsv'"
                shell(
                    "cat {tmpdir}/output/short_summary.specific.{self.lineage}.output.txt > {self.snakemake.output.short_txt:q}"
                )
            if self.snakemake.output.get("short_json"):
                assert self.lineage, "parameter 'lineage' is required to output 'short_json'"
                shell(
                    "cat {tmpdir}/output/short_summary.specific.{self.lineage}.output.json > {self.snakemake.output.short_json:q}"
                )
            if self.snakemake.output.get("full_table"):
                assert self.lineage, "parameter 'lineage' is required to output 'full_table'"
                shell(
                    "cat {tmpdir}/output/run_{self.lineage}/full_table.tsv > {self.snakemake.output.full_table:q}"
                )
            if self.snakemake.output.get("miss_list"):
                assert self.lineage, "parameter 'lineage' is required to output 'miss_list'"
                shell(
                    "cat {tmpdir}/output/run_{self.lineage}/missing_busco_list.tsv > {self.snakemake.output.miss_list:q}"
                )
            if self.snakemake.output.get("out_dir"):
                shell("mv {tmpdir}/output {self.snakemake.output.out_dir:q}")
            if self.snakemake.output.get("dataset_dir"):
                shell("mv {dataset_dir} {self.snakemake.output.dataset_dir:q}")


if __name__ == "__main__":
    Wrapper(snakemake)