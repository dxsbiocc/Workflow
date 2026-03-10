# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/21 16:48:27
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase
from snakemake_wrapper_utils.java import get_java_opts


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.java_opts = get_java_opts(self.snakemake)
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "picard RevertSam"
                " {self.java_opts} {self.extra}"
                " --INPUT {self.snakemake.input[0]}"
                " --TMP_DIR {tmpdir}"
                " --OUTPUT {self.snakemake.output[0]}"
                " {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)