# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/21 14:55:40
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
        control = self.snakemake.input.get("control", None)
        self.control_command = "-i " + control if control else ""

    def run(self):
        shell(
            "(findPeaks"
            " {self.snakemake.input.tag}"
            " -style {self.snakemake.params.style}"
            " {self.extra}"
            " {self.control_command}"
            " -o {self.snakemake.output})"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)