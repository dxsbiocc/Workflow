# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/09 19:15:45
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from snakemake.shell import shell
from tempfile import TemporaryDirectory
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        pass

    def run(self):
        with TemporaryDirectory() as tempdir:
            shell(
                "sambamba sort"
                " {self.snakemake.params}"
                " --nthreads {self.snakemake.threads}"
                " --tmpdir {tempdir}"
                " --out {self.snakemake.output[0]}"
                " {self.snakemake.input[0]}"
                " {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)