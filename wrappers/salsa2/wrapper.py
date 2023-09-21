# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/26 17:22:02
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
        enzyme = self.snakemake.params.get("enzyme", "")
        self.enzyme = f"--enzyme {enzyme}" if enzyme else ""

        gfa = self.snakemake.input.get("gfa", "")
        self.gfa = f"--gfa {gfa}" if gfa else ""
            

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "run_pipeline.py"
                " --assembly {self.snakemake.input.fas}"
                " --length {self.snakemake.input.fai}"
                " --bed {self.snakemake.input.bed}"
                " {self.enzyme}"
                " {self.gfa}"
                " {self.extra}"
                " --output {tmpdir}"
                " {self.log}"
            )

            if self.snakemake.output.get("agp"):
                shell("cat {tmpdir}/scaffolds_FINAL.agp > {self.snakemake.output.agp}")
            if self.snakemake.output.get("fas"):
                shell("cat {tmpdir}/scaffolds_FINAL.fasta > {self.snakemake.output.fas}")


if __name__ == '__main__':
    Wrapper(snakemake)