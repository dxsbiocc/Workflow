# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/26 10:52:23
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
        self.tree = self.snakemake.input.tree
        self.annotation = self.snakemake.input.annotation

        self.output = self.snakemake.output

    def run(self):
        shell(
            "graphlan_annotate.py"
            " --annot {self.annotation}"
            " {self.tree}"
            " {self.output}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)
