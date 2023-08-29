# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/26 17:37:22
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        if self.snakemake.input[0].endswith(".gz"):
            self.pipe = "gunzip -c"
        elif self.snakemake.input[0].endswith(".bz2"):
            self.pipe = "bunzip2 -c"
        else:
            self.pipe = "cat"

    def run(self):
        shell(
            "({self.pipe}"
            " {self.snakemake.input.bedgraph} |"
            " PretextGraph"
            " -i {self.snakemake.input.map_pt}"
            " -n {self.snakemake.params.graph_name}"
            " {self.extra}"
            " -o {self.snakemake.output}"
            ") {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)