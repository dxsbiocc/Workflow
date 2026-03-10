# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/12/09 16:35:32
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from os.path import dirname
from snakemake.shell import shell
from tempfile import TemporaryDirectory

from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):
    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.decoys = self.snakemake.input.get("decoys", "")
        if self.decoys:
            self.decoys = f"--decoys {self.decoys}"

        self.output = self.snakemake.output
        if isinstance(self.output, list) and len(self.output) > 1:
            self.output = dirname(self.snakemake.output[0])

    def run(self):
        with TemporaryDirectory() as tempdir:
            shell(
                "salmon index "
                "--transcripts {self.snakemake.input.sequences} "
                "--index {self.output} "
                "--threads {self.snakemake.threads} "
                "--tmpdir {tempdir} "
                "{self.decoys} "
                "{self.extra} "
                "{self.log}"
            )


if __name__ == "__main__":
    Wrapper(snakemake)