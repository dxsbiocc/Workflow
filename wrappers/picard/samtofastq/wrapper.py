# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 10:52:44
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.java_opts = get_java_opts(self.snakemake)
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)

        fastq1 = self.snakemake.output.fastq1
        fastq2 = self.snakemake.output.get("fastq2", None)
        fastq_unpaired = self.snakemake.output.get("unpaired_fastq", None)

        if not isinstance(fastq1, str):
            raise ValueError("fastq1 needs to be provided")

        self.output = f"--FASTQ {fastq1}"
        if isinstance(fastq2, str):
            self.output += f" --SECOND_END_FASTQ {fastq2}"

        if isinstance(fastq_unpaired, str):
            if not isinstance(fastq2, str):
                raise ValueError("fastq2 is required if fastq_unpaired is set")
            self.output += f" --UNPAIRED_FASTQ {fastq_unpaired}"

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "picard SamToFastq"
                " {self.java_opts} {self.extra}"
                " --INPUT {self.snakemake.input[0]}"
                " --TMP_DIR {tmpdir}"
                " {self.output}"
                " {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)