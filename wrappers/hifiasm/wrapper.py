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
        self.seq = self.snakemake.input.get("seq")

        self.options = ""
        hic1 = self.snakemake.input.get("hic1", None)
        if hic1:
            if isinstance(hic1, list):
                hic1 = " ".join(hic1)
            self.options += " --h1 {}".format(hic1)

        hic2 = self.snakemake.input.get("hic2", None)
        if hic2:
            if isinstance(hic2, list):
                hic2 = " ".join(hic2)
            self.options += " --h2 {}".format(hic2)

        yak = self.snakemake.input.get("yak", None)
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
            " {self.seq}"
            " {self.log}"
        )
        print()
        pattern = re.compile(r'\.bp|\.hic|\.asm')
        for filename in glob.glob(f"{self.output}*.p_ctg.gfa"):
            new_name = pattern.sub('', filename)
            os.rename(filename, new_name)


if __name__ == '__main__':
    Wrapper(snakemake)