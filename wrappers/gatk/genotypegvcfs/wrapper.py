# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/08 10:38:15
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

        intervals = self.snakemake.input.get("intervals", "")
        if not intervals:
            intervals = self.snakemake.params.get("intervals", "")
        self.intervals = "--intervals {}".format(intervals) if intervals else ""

        dbsnp = self.snakemake.input.get("known", "")
        self.dbsnp = "--dbsnp {}".format(dbsnp) if dbsnp else ""

        # Allow for either an input gvcf or GenomicsDB
        gvcf = self.snakemake.input.get("gvcf", "")
        genomicsdb = self.snakemake.input.get("genomicsdb", "")
        if gvcf:
            if genomicsdb:
                raise Exception("Only input.gvcf or input.genomicsdb expected, got both.")
            self.input_string = gvcf
        else:
            if genomicsdb:
                self.input_string = "gendb://{}".format(genomicsdb)
            else:
                raise Exception("Expected input.gvcf or input.genomicsdb.")

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "gatk --java-options '{self.java_opts}' GenotypeGVCFs"
                " --variant {self.input_string}"
                " --reference {self.snakemake.input.ref}"
                " {self.dbsnp}"
                " {self.intervals}"
                " {self.extra}"
                " --tmp-dir {tmpdir}"
                " --output {self.snakemake.output.gvcf}"
                " {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)