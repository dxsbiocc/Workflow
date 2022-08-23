# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 16:40:05
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
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)

        convert_out = self.snakemake.params.get("convert_out", "raw")
        pipes = ""
        if convert_out == "raw":
            pipes = ""
        elif convert_out == "PicardCollectRnaSeqMetrics":
            pipes += " | csvcut -t -c 12,1-10 | csvformat -T"
        else:
            raise ValueError(
                f"Unsupported conversion mode {convert_out}. Please check wrapper documentation."
            )
        self.pipes = pipes

    def run(self):
        shell(
            "(gtfToGenePred"
            " {self.extra}"
            " {self.snakemake.input} /dev/stdout"
            " {self.pipes} > {self.snakemake.output}) {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)