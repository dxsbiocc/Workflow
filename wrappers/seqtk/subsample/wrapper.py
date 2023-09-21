# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/08 16:04:30
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
        self.pairs = True if len(self.snakemake.input) == 2 else False

    def run(self):
        if self.pairs:
            shell(
                "(seqtk sample"
                " -s {self.snakemake.params.seed}"
                " {self.snakemake.input.R1}"
                " {self.snakemake.params.n}"
                " | pigz -9 -p {self.snakemake.threads}"
                " > {self.snakemake.output.R1}"
                " && seqtk sample"
                " -s {self.snakemake.params.seed}"
                " {self.snakemake.input.R2}"
                " {self.snakemake.params.n}"
                " | pigz -9 -p {self.snakemake.threads}"
                " > {self.snakemake.output.R2}"
                ") {self.log}"
            )
        else:
            shell(
                "(seqtk sample"
                " -s {self.snakemake.params.seed}"
                " {self.snakemake.input}"
                " {self.snakemake.params.n}"
                " | pigz -9 -p {self.snakemake.threads}"
                " > {self.snakemake.output}"
                ") {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)