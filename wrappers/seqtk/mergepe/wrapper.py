# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/08 15:56:43
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
        self.compress_lvl = int(self.snakemake.params.get("compress_lvl", 6))

    def run(self):
        shell(
            "(seqtk mergepe {self.snakemake.input} |"
            " pigz -{self.compress_lvl} -c -p {self.snakemake.threads})"
            " > {self.snakemake.output}"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)