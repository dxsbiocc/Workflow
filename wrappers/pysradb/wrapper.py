# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2024/03/28 09:49:46
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import os
import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):
    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.command = self.snakemake.params.get('command')
        assert self.command, "Please specify command"

        self.query = self.snakemake.params.get('query')
        assert self.query, "Please specify query"

    def run(self):
        shell(
            "(pysradb {self.command}"
            " {self.query}"
            " {self.extra}"
            " --saveto {self.snakemake.output}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)