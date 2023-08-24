# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 17:15:27
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
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)
        
        if self.snakemake.output[0].endswith(".gz"):
            self.gzip = True
        self.output = [out.strip(".gz") for out in self.snakemake.output]

        self.options = f" -fq {self.output[0]}"
        if len(self.output) == 2:
            self.options += f" -fq2 {self.output[1]}"


    def run(self):
        shell("(bamToFastq"
        " {self.extra}"
        " -i {self.snakemake.input[0]}"
        " {self.options}"
        ") {self.log}")

        if self.gzip:
            shell("gzip {self.output[0]}")
            if len(self.snakemake.output) == 2:
                shell("gzip {self.output[1]}")


if __name__ == '__main__':
    Wrapper(snakemake)