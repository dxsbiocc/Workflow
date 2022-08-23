# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/18 20:48:44
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import os
import pathlib
from tempfile import TemporaryDirectory

from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        assert isinstance(self.snakemake.input, str), "Only supported 'str'!"

    def run(self):
        # Run fastqc, since there can be race conditions if multiple jobs
        # use the same fastqc dir, we create a temp dir.
        with TemporaryDirectory() as tempdir:
            shell(
                "fastqc {self.snakemake.params} -t {self.snakemake.threads} "
                "--outdir {tempdir:q} {self.snakemake.input[0]:q}"
                " {self.log}"
            )

            # Move outputs into proper position.
            output_base = pathlib.Path(self.snakemake.input[0]).stem
            html_path = os.path.join(tempdir, f"{output_base}_fastqc.html")
            zip_path = os.path.join(tempdir, f"{output_base}_fastqc.zip")

            if self.snakemake.output.html != html_path:
                shell("mv {html_path:q} {self.snakemake.output.html:q}")

            if self.snakemake.output.zip != zip_path:
                shell("mv {zip_path:q} {self.snakemake.output.zip:q}")


if __name__ == '__main__':
    Wrapper(snakemake)