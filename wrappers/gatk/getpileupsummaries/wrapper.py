# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/08 10:38:32
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
        with tempfile.TemporaryDirectory() as tempdir:
            shell(
                "gatk GetPileupSummaries "
                "--java-options '{self.java_opts}' "
                "--input {self.snakemake.input.bam} "
                "--intervals {self.snakemake.input.intervals} "
                "--variant {self.snakemake.input.variants} "
                "--output {self.snakemake.output[0]} "
                "--tmp-dir {tempdir} "
                "{self.extra} "
                "{self.log} "
            )


if __name__ == '__main__':
    Wrapper(snakemake)