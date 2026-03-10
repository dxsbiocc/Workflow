# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/08 09:51:40
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

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "gatk --java-options '{self.java_opts}' ApplyVQSR"
                " --variant {self.snakemake.input.vcf}"
                " --recal-file {self.snakemake.input.recal}"
                " --reference {self.snakemake.input.ref}"
                " --tranches-file {self.snakemake.input.tranches}"
                " --mode {self.snakemake.params.mode}"
                " {self.extra}"
                " --tmp-dir {tmpdir}"
                " --output {self.snakemake.output.vcf}"
                " {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)