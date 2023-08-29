# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/26 14:35:47
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import re
import os
import glob

from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.db = os.path.abspath(self.snakemake.input.get("db"))

        self.script = os.path.abspath(self.snakemake.input.get("script"))

        prefix = os.path.abspath(self.snakemake.input.get("seq", 'all'))
        self.prefix =os.path.splitext(prefix)[0]

        self.sample = os.path.basename(self.prefix)
    
        self.min_length = self.snakemake.params.get("min_length", 44)
        self.min_percentage = self.snakemake.params.get("min_percentage", 97)


    def run(self):
        shell(
            "({self.script}"
            " -p {self.prefix} "
            " -d {self.db}"
            " -l {self.min_length}"
            " -m {self.min_percentage}"
            " -t {self.snakemake.threads}"
            " -o {self.snakemake.output}"
            ") {self.log}"
        )
        shell("mv {self.snakemake.output}/{self.sample}.stats {self.snakemake.output}/{self.sample}.hifi.stats")


if __name__ == '__main__':
    Wrapper(snakemake)