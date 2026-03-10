# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2024/06/22 14:12:06
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
        fq1 = self.snakemake.input.get("fq1")                        # Path to the FASTQ R1 file(s)
        fq2 = self.snakemake.input.get("fq2")                        # Path to the FASTQ R2 file(s)
        self.fastq = ' '.join([f'{r1} {r2}' for r1, r2 in zip(fq1, fq2)])
        # required arguments
        index = self.snakemake.input.get("index")                         # Path to the kallisto index
        self.index = os.path.join(index, 'kb_ref.idx')
        self.t2g = os.path.join(index, 't2g.txt')                         # Path to transcript-to-gene mapping

        self.tech = self.snakemake.params.get("tech", "10x")              # Sequencing technology

        self.mode = self.snakemake.params.get("mode", "standard")
        if self.mode not in ["standard", "nac", "kite", "custom"]:
            raise ValueError("Invalid mode: {}".format(self.mode))
        
        # required arguments for `nac` workflow
        self.cdna, self.intron = "", ""
        if self.mode == 'nac':
            self.cdna = f"-c1 {os.path.join(index, 'cdna_t2c.txt')}"      # Path to cDNA transcripts-to-capture
            self.intron = f"-c2 {os.path.join(index, 'intron_t2c.txt')}"  # Path to intron transcripts-to-captured

        white = self.snakemake.input.get("whitelist")                     # Path to the whitelist file
        self.whitelist = f"-w {white}" if white else ""

    def run(self):
        shell(
            "kb count"
            " -t {self.snakemake.threads}"
            " -i {self.index}"
            " -g {self.t2g}"
            " {self.whitelist}"
            " {self.cdna} {self.intron}"
            " -x {self.tech}"
            " --workflow {self.mode}"
            " {self.extra}"
            " -o {self.snakemake.output}"
            " {self.fastq}"
            " {self.log}"
        )
    

if __name__ == "__main__":
    Wrapper(snakemake)