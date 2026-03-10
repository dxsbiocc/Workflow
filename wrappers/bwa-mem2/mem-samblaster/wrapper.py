# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/29 20:40:05
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
        # Extract arguments.
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)
        self.sort_extra = self.snakemake.params.get("sort_extra", "")
        self.samblaster_extra = self.snakemake.params.get("samblaster_extra", "")

        index = self.snakemake.input.get("index", "")
        if isinstance(index, str):
            self.index = path.splitext(self.snakemake.input.idx)[0]
        else:
            self.index = path.splitext(self.snakemake.input.idx[0])[0]


        # Check inputs/arguments.
        if not isinstance(self.snakemake.input.reads, str) and len(self.snakemake.input.reads) not in {
            1,
            2,
        }:
            raise ValueError("input must have 1 (single-end) or 2 (paired-end) elements")

    def run(self):
        shell(
            "(bwa-mem2 mem"
            " -t {self.snakemake.threads}"
            " {self.extra}"
            " {self.index}"
            " {self.snakemake.input.reads}"
            " | samblaster"
            " {self.samblaster_extra}"
            " | sambamba view -S -f bam /dev/stdin"
            " -t {self.snakemake.threads}"
            " | sambamba sort /dev/stdin"
            " -t {self.snakemake.threads}"
            " -o {self.snakemake.output.bam}"
            " {self.sort_extra}"
            ") {self.log}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)