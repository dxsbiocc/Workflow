# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 09:29:47
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import os
import sys
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        percentage = self.snakemake.params.get('percentage')
        number = self.snakemake.params.get('number')

        assert bool(percentage) ^ bool(number), \
            "Parameters `-p` and `-n` can not setted at the same time"
        
        if percentage:
            self.extra += f' -p {percentage}'
        if number:
            self.extra += f' -n {number}'


    def run(self):
        shell(
            "(macs2 randsample "
            " -i {self.snakemake.input}"
            " -o {self.snakemake.output}"
            " {self.extra}"
            ") {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)