# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/08 13:11:23
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

        vcf = self.snakemake.input.vcf
        if isinstance(vcf, str):
            self.vcf = "--eval  {}".format(vcf)
        else:
            self.vcf = list(map("--eval {}".format, vcf))

        bam = self.snakemake.input.get("bam", "")
        if bam:
            if isinstance(bam, str):
                self.bam = "--input  {}".format(bam)
            else:
                self.bam = list(map("--input {}".format, bam))
        else:
            self.bam = ""

        ref = self.snakemake.input.get("ref", "")
        self.ref = "--reference " + ref if ref else ""
            
        ref_dict = self.snakemake.input.get("dict", "")
        self.ref_dict = "--sequence-dictionary " + ref_dict if ref_dict else ""
            
        known = self.snakemake.input.get("known", "")
        self.known = "--dbsnp " + known if known else ""
            
        comp = self.snakemake.input.get("comp", "")
        if comp:
            if isinstance(comp, str):
                self.comp = "--comparison  {}".format(comp)
            else:
                self.comp = list(map("--comparison {}".format, comp))
        else:
            self.comp = ""

        ped = self.snakemake.input.get("ped", "")
        self.ped = "--pedigree " + ped if ped else ""
            
    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "gatk --java-options '{self.java_opts}' VariantEval"
                " {self.vcf}"
                " {self.bam}"
                " {self.ref}"
                " {self.ref_dict}"
                " {self.known}"
                " {self.ped}"
                " {self.comp}"
                " {self.extra}"
                " --tmp-dir {tmpdir}"
                " --output {self.snakemake.output[0]}"
                " {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)