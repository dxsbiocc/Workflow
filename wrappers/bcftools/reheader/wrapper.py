# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/21 16:18:22
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import pathlib
import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase
from snakemake_wrapper_utils.bcftools import get_bcftools_opts


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)
        # bcftools options
        self.bcftools_opts = get_bcftools_opts(
            self.snakemake, parse_ref=False, parse_memory=False
        )
        self.view_extra = self.snakemake.params.get("view_extra", "")

        ## Extract arguments: header
        header = self.snakemake.input.get("header", None)
        self.header = f"-h {header}" if header else ""
        # samples
        samples = self.snakemake.input.get("samples", "")
        self.samples = f"-s {samples}" if samples else ""

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_prefix = pathlib.Path(tmpdir) / "bcftools_reheader."

            shell(
                "(bcftools reheader"
                " --threads {self.snakemake.threads}"
                " {self.header}"
                " {self.samples}"
                " {self.extra}"
                " --temp-prefix {tmp_prefix}"
                " {self.snakemake.input[0]}"
                "| bcftools view"
                " {self.bcftools_opts}"
                " {self.view_extra}"
                ") {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)