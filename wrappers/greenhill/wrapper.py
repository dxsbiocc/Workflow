# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/09/26 18:51:52
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
        aware = self.snakemake.input.get("aware", "")
        cph = self.snakemake.input.get("cph", "")
        assert not (aware and cph), "'aware' and 'cph' can not be used at the same time"
        if aware:
            self.extra += f" -c {aware}"
            self.extra += ' -b {}'.format(' '.join(self.snakemake.input.get("bubble", [])))
        if cph:
            self.extra += f" -cph {cph}"
        
        long = self.snakemake.input.get("long", "")
        if long:
            if isinstance(long, list):
                long = ' '.join(long)
            self.extra += f" -p {long}"

        tmp = self.snakemake.params.get("tmp", "")
        if tmp:
            self.extra += f" -tmp {tmp}"

        hic = self.snakemake.params.get("hic", "paired")
        assert hic in ["paired", "single"], "'hic' must be 'paired' or 'single'"
        if hic == 'paired':
            self.hic = 'HIC'
        else:
            self.hic = 'hic'

        output = os.path.abspath(self.snakemake.output[0])
        prefix = self.snakemake.params.get("prefix", "out")
        self.output = os.path.join(output, prefix)

    def run(self):
        shell(
            "greenhill"
            " {self.extra}"
            " -t {self.snakemake.threads}"
            " {self.hic} {self.snakemake.input.hic}"
            " -o {self.output}"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)