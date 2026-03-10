# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2024/01/11 20:01:59
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file
# ============================================================


import os
from snakemake import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.input = self.snakemake.input
        if isinstance(self.input, str):
            self.input = [self.input]
        self.output = self.snakemake.output
        # parameters
        self.sample_list = self.snakemake.params.get("sample_list")
        self.metrics = self.snakemake.params.get("metrics")

        if isinstance(self.sample_list, str):
            self.sample_list = [self.sample_list]
        if isinstance(self.metrics, str):
            self.metrics = [self.metrics]

        assert len(self.sample_list) == len(
            self.input
        ), "The number of samples must be equal to the number of input files"

    def run(self):
        shell(f"echo -e 'SampleID\t{'\t'.join(self.metrics)}' > {self.output}")
        for i, infile in enumerate(self.input):
            shell(f"echo -e -n '{self.sample_list[i]}\t' >> {self.output}")
            for metric in self.metrics:
                shell(
                    f"alpha_diversity.py -f {infile} -a {metric}"
                    " | cut -f 2 -d ':' | tr '\n' '\t'"
                    f" >> {self.output}"
                )
            shell("echo '' >> {self.output}")


if __name__ == "__main__":
    Wrapper(snakemake)
