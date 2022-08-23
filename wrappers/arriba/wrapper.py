# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/23 09:38:54
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
        discarded_fusions = self.snakemake.output.get("discarded", "")
        self.discarded_cmd = f"-O {discarded_fusions}" if discarded_fusions else ""
        # blacklist
        blacklist = self.snakemake.params.get("blacklist")
        self.blacklist_cmd = f"-b {blacklist}" if blacklist else ""

        known_fusions = self.snakemake.params.get("known_fusions")
        self.known_cmd = f"-k {known_fusions}" if known_fusions else ""

        sv_file = self.snakemake.params.get("sv_file")
        self.sv_cmd = f"-d {sv_file}" if sv_file else ""

    def run(self):
        shell(
            "arriba "
            "-x {self.snakemake.input.bam} "
            "-a {self.snakemake.input.genome} "
            "-g {self.snakemake.input.annotation} "
            "{self.blacklist_cmd} "
            "{self.known_cmd} "
            "{self.sv_cmd} "
            "-o {self.snakemake.output.fusions} "
            "{self.discarded_cmd} "
            "{self.extra} "
            "{self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)