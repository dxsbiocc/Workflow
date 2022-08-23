# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 11:43:13
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import pathlib
import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.samtools import get_samtools_opts
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.samtools_opts = get_samtools_opts(self.snakemake)

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_prefix = pathlib.Path(tmpdir) / "samtools_fastq.sort_"
            shell(
                "samtools sort {self.samtools_opts} {self.extra} -T {tmp_prefix} {self.snakemake.input[0]} {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)