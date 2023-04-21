# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/18 20:07:45
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import re
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.adapters = self.snakemake.params.get("adapters", "")
        # Assert input files
        n = len(self.snakemake.input.reads)
        if n == 1:
            self.reads = "--in1 {}".format(self.snakemake.input.reads)
        elif n == 2:
            self.reads = "--in1 {} --in2 {}".format(*self.snakemake.input.reads)
        else:
            raise ValueError("input->sample must have 1 (single-end) or 2 (paired-end) elements.")

        # Output files
        trimmed_paths = self.snakemake.output.get("trimmed", None)
        if trimmed_paths:
            if n == 1:
                trimmed = "--out1 {}".format(self.snakemake.output.trimmed)
            else:
                trimmed = "--out1 {} --out2 {}".format(*self.snakemake.output.trimmed)
                # Output unpaired files
                unpaired = self.snakemake.output.get("unpaired", None)
                if unpaired:
                    trimmed += f" --unpaired1 {unpaired} --unpaired2 {unpaired}"
                else:
                    unpaired1 = self.snakemake.output.get("unpaired1", None)
                    if unpaired1:
                        trimmed += f" --unpaired1 {unpaired1}"
                    unpaired2 = self.snakemake.output.get("unpaired2", None)
                    if unpaired2:
                        trimmed += f" --unpaired2 {unpaired2}"

                # Output merged PE reads
                merged = self.snakemake.output.get("merged", None)
                if merged:
                    if not re.search(r"--merge\b", self.extra):
                        raise ValueError(
                            "output.merged specified but '--merge' option missing from params.extra"
                        )
                    trimmed += f" --merged_out {merged}"
        else:
            trimmed = ""
        
        # Output failed reads
        failed = self.snakemake.output.get("failed", None)
        if failed:
            trimmed += f" --failed_out {failed}"
        self.trimmed = trimmed
        # Stats
        self.html = "--html {}".format(self.snakemake.output.html)
        self.json = "--json {}".format(self.snakemake.output.json)

    def run(self):
        shell(
            "(fastp --thread {self.snakemake.threads} "
            "{self.extra} "
            "{self.adapters} "
            "{self.reads} "
            "{self.trimmed} "
            "{self.json} "
            "{self.html} ) {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)