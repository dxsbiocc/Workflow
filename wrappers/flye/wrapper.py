# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/09/02 09:27:09
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import os
import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):
    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.mode = self.snakemake.params.get("mode", "")
        assert self.mode in [
            "--pacbio-raw", 
            "--pacbio-corr", 
            "--pacbio-hifi", 
            "--nano-raw", 
            "--nano-corr", 
            "--nano-hq"
        ], "Unrecognised mode to run Flye. Options: --pacbio-raw, --pacbio-corr, --pacbio-hifi, --nano-raw, --nano-corr, --nano-hq."
        
        self.assembly = self.snakemake.output.assembly
        self.assembly_gfa = self.snakemake.output.assembly_gfa
        self.assembly_gv = self.snakemake.output.assembly_gv

        for _, value in self.snakemake.output.items():
            dirname = os.path.dirname(value)
            if not os.path.exists(dirname):
                os.makedirs(dirname)

    
    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "(flye {self.mode}"
                " {self.snakemake.input}"
                " --out-dir {tmpdir}"
                " --threads {self.snakemake.threads}"
                " {self.extra}"
                ") {self.log}"
            )
            
            if self.assembly:
                if self.assembly.endswith(".gz"):
                    shell("gzip -c {tempdir}/assembly.fasta > {self.assembly}")
                else:
                    shell("mv {tempdir}/assembly.fasta {self.assembly}")

            if self.assembly_gfa:
                if self.assembly_gfa.endswith(".gz"):
                    shell("gzip -c {tempdir}/assembly_graph.gfa > {self.assembly_gfa}")
                else:
                    shell("mv {tempdir}/assembly_graph.gfa {self.assembly_gfa}")

            if self.assembly_gv:
                if self.assembly_gv.endswith(".gz"):
                    shell("gzip -c {tempdir}/assembly_graph.gv > {self.assembly_gv}")
                else:
                    shell("mv {tempdir}/assembly_graph.gv {self.assembly_gv}")

            # shell('mv {tempdir}/flye.log {self.snakemake.log}')
            shell('mv {tempdir}/assembly_info.txt {self.snakemake.output.assembly_info}')
            shell('mv {tempdir}/params.json {self.snakemake.output.params}')


if __name__ == "__main__":
    wrapper = Wrapper(snakemake)