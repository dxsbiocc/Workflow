# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/08 12:33:09
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.java_opts = get_java_opts(self.snakemake)
        if self.snakemake.output.get("bam", None):
            self.bam_output = "--bam-output " + self.snakemake.output.bam
        else:
            self.bam_output = ""

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "gatk --java-options '{self.java_opts}' Mutect2"  # Tool and its subprocess
                " --native-pair-hmm-threads {self.snakemake.threads}"
                " --input {self.snakemake.input.map}"  # Path to input mapping file
                " --reference {self.snakemake.input.fasta}"  # Path to reference fasta file
                " {self.extra}"  # Extra parameters
                " --tmp-dir {tmpdir}"
                " --output {self.snakemake.output.vcf}"  # Path to output vcf file
                " {self.bam_output}"  # Path to output bam file, optional
                " {self.log}"  # Logging behaviour
            )


if __name__ == '__main__':
    Wrapper(snakemake)