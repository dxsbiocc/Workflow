# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/08 09:55:16
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

        self.known = self.snakemake.input.get("known", "")
        if self.known:
            if isinstance(self.known, str):
                self.known = [self.known]
            self.known = list(map("--known-sites {}".format, self.known))

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "gatk --java-options '{self.java_opts}' BaseRecalibrator"
                " --input {self.snakemake.input.bam}"
                " --reference {self.snakemake.input.ref}"
                " {self.known}"
                " {self.extra}"
                " --tmp-dir {tmpdir}"
                " --output {self.snakemake.output.recal_table}"
                " {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)