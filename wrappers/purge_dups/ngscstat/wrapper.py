# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/28 10:52:20
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        pass

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell("ngscstat {self.extra} -O {tmpdir} {self.snakemake.input} {self.log}")

            if self.snakemake.output.get("cov"):
                shell("cat {tmpdir}/TX.base.cov > {self.snakemake.output.cov}")
            if self.snakemake.output.get("stat"):
                shell("cat {tmpdir}/TX.stat > {self.snakemake.output.stat}")


if __name__ == '__main__':
    Wrapper(snakemake)