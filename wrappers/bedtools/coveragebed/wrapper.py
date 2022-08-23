# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 19:49:49
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
        shell.executable("bash")

        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)

        self.input_a = self.snakemake.input.a
        self.input_b = self.snakemake.input.b

        self.output_file = self.snakemake.output[0]

        if not isinstance(self.output_file, str) and len(self.snakemake.output) != 1:
            raise ValueError("Output should be one file: " + str(self.output_file) + "!")

        shell(
            "coverageBed"
            " -a {self.input_a}"
            " -b {self.input_b}"
            " {self.extra}"
            " > {self.output_file}"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)