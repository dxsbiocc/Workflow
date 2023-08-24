# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/17 20:37:26
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
        self.fastq = self.snakemake.input.get("fastq")

        self.options = ""
        hic = self.snakemake.params.get("hic", None)
        if hic and len(hic) == 2:
            self.options += " --h1 {hic[0]} --h2 {hic[1]}"

        yak = self.snakemake.params.get("yak", None)
        if yak and len(yak) == 2:
            self.options += " -1 {yak[0]} -2 {yak[1]}"

        ultralong = self.snakemake.params.get("ultralong")
        if ultralong:
            self.options += " --ul {ultralong}"

        if self.snakemake.params.get("primary", False):
            self.options += " --primary"

        self.output = os.path.abspath(os.path.commonprefix(self.snakemake.output).strip('.'))

    def run(self):
        shell(
            "hifiasm -t {self.snakemake.threads} "
            " -o {self.output}"
            " {self.options}"
            " {self.extra}"
            " {self.fastq}"
            " {self.log}"
        )
        pattern = re.compile(r'\.bp|\.hic|\.asm')
        for filename in glob.glob(f"{self.output}*.p_ctg.gfa"):
            new_name = pattern.sub('', filename)
            os.rename(filename, new_name)


if __name__ == '__main__':
    Wrapper(snakemake)