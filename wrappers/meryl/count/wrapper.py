# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/28 14:25:19
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase
from snakemake_wrapper_utils.snakemake import get_mem


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.command = self.snakemake.params.get("command", "count")
        assert self.command in [
            "count",
            "count-forward",
            "count-reverse",
        ], "invalid command specified."


        self.mem_gb = get_mem(self.snakemake, out_unit="GiB")

    def run(self):
        shell(
            "meryl"
            " {self.command}"
            " threads={self.snakemake.threads}"
            " memory={self.mem_gb}"
            " {self.extra}"
            " {self.snakemake.input}"
            " output {self.snakemake.output}"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)