# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 14:09:10
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
        # Samtools takes additional threads through its option -@
        # One thread for samtools merge
        # Other threads are *additional* threads passed to the '-@' argument
        self.threads = "" if self.snakemake.threads <= 1 else \
            " -@ {} ".format(self.snakemake.threads - 1)

    def run(self):
        shell(
            "samtools index"
            " {self.threads} {self.extra} "
            " {self.snakemake.input[0]}"
            " {self.snakemake.output[0]}"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)