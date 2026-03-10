# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/28 14:42:45
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
        self.command = self.snakemake.params.get("command", "union")
        assert self.command in [
            "union",
            "union-min",
            "union-max",
            "union-sum",
            "intersect",
            "intersect-min",
            "intersect-max",
            "intersect-sum",
            "subtract",
            "difference",
            "symmetric-difference",
        ], "invalid command specified."

    def run(self):
        shell(
            "meryl {self.command}"
            " {self.snakemake.input}"
            " output {self.snakemake.output}"
            " {self.log}")


if __name__ == '__main__':
    Wrapper(snakemake)