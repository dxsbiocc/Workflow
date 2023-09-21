# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/26 17:42:42
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import tempfile
from pathlib import Path
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.sequences = "=" + ", =".join(self.snakemake.output.keys())


        format = Path(self.snakemake.output.full).suffix.lstrip(".")
        self.format = "jpeg" if format == "jpg" else format

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "PretextSnapshot"
                " --map {self.snakemake.input[0]}"
                " {self.extra}"
                " --sequences {self.sequences:q}"
                " --format {self.format}"
                " --folder {tmpdir}"
                " --prefix out_"
                " {self.log}"
            )

            if self.snakemake.output.get("full"):
                shell("mv {tmpdir}/out_FullMap.{self.format} {self.snakemake.output.full}")
            if self.snakemake.output.get("all"):
                Path(self.snakemake.output.all).mkdir(parents=True, exist_ok=True)
                shell("mv {tmpdir}/out_*.{self.format} {self.snakemake.output.all}/")


if __name__ == '__main__':
    Wrapper(snakemake)