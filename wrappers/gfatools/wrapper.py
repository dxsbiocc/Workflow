# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/26 10:52:23
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.command = self.snakemake.params.get("command", "view")
        assert self.command in [
            "view",
            "stat",
            "gfa2fa",
            "gfa2bed",
            "blacklist",
            "bubble",
            "asm",
            "sql",
            "version",
        ], "invalid command specified."

    def run(self):
        shell(
            "gfatools"
            " {self.command}"
            " {self.extra}"
            " {self.snakemake.input}"
            " > {self.snakemake.output}"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)