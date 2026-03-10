# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/12/09 16:39:35
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file
# ============================================================


from os.path import dirname, isfile, isdir
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):
    def __init__(self, snakemake):
        super().__init__(snakemake)

    def parser(self):
        quants = self.snakemake.input.get("quants", "")
        if isinstance(quants, str):
            quants = [quants]
        self.quants = []
        for quant in quants:
            if isfile(quant):
                self.quants.append(dirname(quant))
            elif isdir(quant):
                self.quants.extend(quant)
            else:
                raise ValueError(f"The input {quant} is not a file or directory!")

        self.output = self.snakemake.output

    def run(self):
        shell(
            "salmon quantmerge"
            " --quants {self.quants} "
            " {self.extra}"
            " -o {self.output}"
            " {self.log}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)
