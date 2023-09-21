# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 21:39:10
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
        genome = self.snakemake.input.get("genome", None)
        faidx = self.snakemake.input.get("faidx", None)

        if bool(genome) ^ bool(faidx):
            if genome:
                self.extra += " -g {}".format(genome)
            else:
                self.extra += " -faidx {}".format(faidx)
        else:
            assert not genome, "sort file either by the variable genome or faidx, not both."

    def run(self):
        shell(
            "(bedtools sort"
            " {self.extra}"
            " -i {self.snakemake.input.infile}"
            " > {self.snakemake.output[0]})"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)