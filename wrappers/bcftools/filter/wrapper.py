# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/21 15:34:55
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase
from snakemake_wrapper_utils.bcftools import get_bcftools_opts


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.bcftools_opts = get_bcftools_opts(self.snakemake, parse_ref=False, parse_memory=False)
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)
        self.filter = self.snakemake.params.get("filter", "")


        if isinstance(self.snakemake.output, list) and len(self.snakemake.output) > 1:
            raise Exception("Only one output file expected, got: " + str(len(self.snakemake.output)))

    def run(self):
        shell(
            "bcftools filter"
            " {self.bcftools_opts}"
            " {self.filter}"
            " {self.extra}"
            " {self.snakemake.input[0]}"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)