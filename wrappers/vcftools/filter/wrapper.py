# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/20 21:01:14
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.input_flag = "--vcf"
        if self.snakemake.input[0].endswith(".gz"):
            self.input_flag = "--gzvcf"

        self.output = f" > {self.snakemake.output[0]}"
        if self.output.endswith(".gz"):
            self.output = " | gzip" + self.output

        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)

    def run(self):
        shell(
            "(vcftools "
            " {self.input_flag}"
            " {self.snakemake.input}"
            " {self.extra}"
            " --recode"
            " --stdout"
            " {self.output}"
            ") {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)