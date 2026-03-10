# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/29 16:08:54
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
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)

        agpfile = self.snakemake.params.get("agpfile", "")
        if agpfile:
            self.extra += f" --agp-to-path {agpfile}"
        include_bed = self.snakemake.params.get("include_bed", "")
        if include_bed:
            self.extra += f" --include-bed {include_bed}"
        exclude_bed = self.snakemake.params.get("exclude_bed", "")
        if exclude_bed:
            self.extra += f" --exclude-bed {exclude_bed}"
        instructions = self.snakemake.params.get("instructions", "")
        if instructions:
            self.extra += f" --swiss-army-knife {instructions}"

    def run(self):
        shell(
            "gfastats"
            " {self.extra}"
            " --threads {self.snakemake.threads}"
            " --out-format {self.snakemake.output.seq}"
            " --input-sequence {self.snakemake.input}"
            " > {self.snakemake.output.summary}"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)