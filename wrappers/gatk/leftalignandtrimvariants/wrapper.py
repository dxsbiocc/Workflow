# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/08 15:24:14
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
            

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "gatk --java-options '{self.java_opts}' LeftAlignAndTrimVariants"
                " --variant {self.snakemake.input.vcf}"
                " --reference {self.snakemake.input.ref}"
                " {self.intervals}"
                " {self.extra}"
                " --tmp-dir {tmpdir}"
                " --output {self.snakemake.output.vcf}"
                " {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)