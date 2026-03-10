# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/08 10:27:34
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
        self.aln = self.snakemake.input.get("aln", "")
        if self.aln:
            self.aln = f"--input {self.aln}"

        self.contamination = self.snakemake.input.get("contemination_table", "")
        if self.contamination:
            self.contamination = f"--contamination-table {self.contamination}"

        self.segmentation = self.snakemake.input.get("segmentation", "")
        if self.segmentation:
            self.segmentation = f"--tumor-segmentation {self.segmentation}"

        self.f1r2 = self.snakemake.input.get("f1r2", "")
        if self.f1r2:
            self.f1r2 = f"--orientation-bias-artifact-priors {self.f1r2}"

        self.intervals = self.snakemake.input.get("bed", "")
        if self.intervals:
            self.intervals = f"--intervals {self.intervals}"

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "gatk --java-options '{self.java_opts}' FilterMutectCalls"
                " --variant {self.snakemake.input.vcf}"
                " --reference {self.snakemake.input.ref}"
                " {self.aln}"  # BAM/SAM/CRAM file containing reads
                " {self.contamination}"  # Tables containing contamination information
                " {self.segmentation}"  # Tumor segments' minor allele fractions
                " {self.f1r2}"  # .tar.gz files containing tables of prior artifact
                " {self.intervals}"  # Genomic intervals over which to operate
                " {self.extra}"
                " --tmp-dir {tmpdir}"
                " --output {self.snakemake.output.vcf}"
                " {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)