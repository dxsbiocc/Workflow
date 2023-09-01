# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/26 17:38:06
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase
from snakemake_wrapper_utils.samtools import infer_out_format


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        input_file = self.snakemake.input[0]
        if input_file.endswith(".gz"):
            self.pipe = "gunzip -c"
        elif input_file.endswith(".bz2"):
            self.pipe = "bunzip2 -c"
        elif infer_out_format(input_file) in ["SAM", "BAM", "CRAM"]:
            self.pipe = "samtools view -h"
        else:
            self.pipe = "cat"

    def run(self):
        shell(
            "({self.pipe}"
            " {self.snakemake.input} | "
            "PretextMap"
            " {self.extra}"
            " -o {self.snakemake.output}"
            ") {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)