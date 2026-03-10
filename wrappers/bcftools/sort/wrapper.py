# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/21 16:24:42
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase
from snakemake_wrapper_utils.bcftools import get_bcftools_opts


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.bcftools_opts = get_bcftools_opts(self.snakemake, parse_ref=False)

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "bcftools sort"
                " {self.bcftools_opts}"
                " {self.extra}"
                " --temp-dir {tmpdir}"
                " {self.snakemake.input[0]}"
                " {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)