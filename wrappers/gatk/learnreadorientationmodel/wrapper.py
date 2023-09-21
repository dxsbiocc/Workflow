# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/08 11:19:08
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

        self.f1r2 = "--input "
        if isinstance(self.snakemake.input["f1r2"], list):
            # Case user provided a list of archives
            self.f1r2 += "--input ".join(self.snakemake.input["f1r2"])
        else:
            # Case user provided a single archive as a string
            self.f1r2 += self.snakemake.input["f1r2"]

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "gatk --java-options '{self.java_opts}' LearnReadOrientationModel"  # Tool and its subprocess
                " {self.f1r2}"  # Path to input mapping file
                " {self.extra}"  # Extra parameters
                " --tmp-dir {tmpdir}"
                " --output {self.snakemake.output[0]}"  # Path to output vcf file
                " {self.log}"  # Logging behaviour
            )


if __name__ == '__main__':
    Wrapper(snakemake)