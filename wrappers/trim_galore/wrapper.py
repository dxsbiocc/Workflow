# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/09 09:04:04
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import os
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        # Check that two input files were supplied
        n = len(self.snakemake.input)
        self.paired = "--paired" if n == 2 else ""

        # Don't run with `--fastqc` flag
        if "--fastqc" in self.snakemake.params.get("extra", ""):
            raise ValueError(
                "The trim_galore Snakemake wrapper cannot "
                "be run with the `--fastqc` flag. Please "
                "remove the flag from extra params. "
                "You can use the fastqc Snakemake wrapper on "
                "the input and output files instead."
            )

        # Check that four output files were supplied
        output = self.snakemake.output.trimmed + self.snakemake.output.report
        m = len(output)
        assert m == 2 * n, "Output must contain 2/4 files. Given: %r." % m

        # Check that all output files are in the same directory
        self.out_dir = os.path.dirname(output[0])
        for file_path in output[1:]:
            assert self.out_dir == os.path.dirname(file_path), (
                "trim_galore can only output files to a single directory."
                " Please indicate only one directory for the output files."
            )

    def run(self):
        shell(
            "(trim_galore"
            " {self.extra}"
            " {self.paired}"
            " -o {self.out_dir}"
            " {self.snakemake.input})"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)