# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/09 16:39:03
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
        include_bed = self.snakemake.input.get("include_bed", "")
        exclude_bed = self.snakemake.input.get("exclude_bed", "")

        if include_bed and exclude_bed:
            raise Exception("supply either include_bed or exclude_bed, not both")

        self.include_bed = f"-j {include_bed}" if include_bed else ""
        self.exclude_bed = f"-J {exclude_bed}" if exclude_bed else ""

        self.output_prefix = self.snakemake.output[0].rsplit("_", 1)[0]

    def run(self):
        shell(
            "pindel "
            "-T {self.snakemake.threads} "
            "{self.extra} "
            "{self.include_bed} "
            "{self.exclude_bed} "
            "-i {self.snakemake.input.config} "
            "-f {self.snakemake.input.ref} "
            "-o {self.output_prefix} {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)