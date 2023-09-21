# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/08 10:38:41
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
        bams = self.snakemake.input.bam
        if isinstance(bams, str):
            bams = [bams]
        self.bams = list(map("--input {}".format, bams))

        intervals = self.snakemake.input.get("intervals", "")
        if not intervals:
            intervals = self.snakemake.params.get("intervals", "")
        self.intervals = "--intervals {}".format(intervals) if intervals else ""

        known = self.snakemake.input.get("known", "")
        self.known = "--dbsnp " + str(known) if known else ""

        vcf_output = self.snakemake.output.get("vcf", "")
        if vcf_output:
            self.output = " --output " + str(vcf_output)

        gvcf_output = self.snakemake.output.get("gvcf", "")
        if gvcf_output:
            self.output = " --emit-ref-confidence GVCF " + " --output " + str(gvcf_output)

        if (vcf_output and gvcf_output) or (not gvcf_output and not vcf_output):
            if vcf_output and gvcf_output:
                raise ValueError(
                    "please set vcf or gvcf as output, not both! It's not supported by gatk"
                )
            else:
                raise ValueError("please set one of vcf or gvcf as output (not both)!")

        bam_output = self.snakemake.output.get("bam", "")
        self.bam_output = " --bam-output " + str(bam_output) if bam_output else ""

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "gatk --java-options '{self.java_opts}' HaplotypeCaller"
                " --native-pair-hmm-threads {self.snakemake.threads}"
                " {self.bams}"
                " --reference {self.snakemake.input.ref}"
                " {self.intervals}"
                " {self.known}"
                " {self.extra}"
                " --tmp-dir {tmpdir}"
                " {self.output}"
                " {self.bam_output}"
                " {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)