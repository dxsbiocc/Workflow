# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/30 20:06:49
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from os import path
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.command = self.snakemake.params.get("command", "")
        assert self.command in [
            "hiclib",
            "pairix",
            "pairs",
            "tabix"
        ], "command must be one of [hiclib, pairix, pairs, tabix]"

        bins = self.snakemake.input.get("bins", "")
        binsize = self.snakemake.params.get("binsize", None)
        if bins.endswith(".bed"):
            self.bins = bins
        elif bins and binsize:
            self.bins = f" {bins}:{binsize}"
        else:
            raise ValueError("bins must be a bed file or binsize must be specified")


    def run(self):
        shell(
            "(cooler cload {self.command}"
            " {self.extra}"
            " {self.bins}"
            " {self.snakemake.input.pairs}"
            " {self.snakemake.output}"
            ") {self.log}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)