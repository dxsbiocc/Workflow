# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/18 15:29:35
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
        reads = self.snakemake.input.get("reads", "")
        if reads:
            if isinstance(reads, list):
                self.extra = f" -c {' '.join(reads)}"
            elif isinstance(reads, str):
                self.extra = f" --single {reads}"
            else:
                raise ValueError(
                    "Input must be a list of fastq files or a single fastq file."
                )
        else:
            r1 = self.snakemake.input.get("r1", "")
            r2 = self.snakemake.input.get("r2", "")
            if r1 and r2:
                self.extra = f" -1 {r1} -2 {r2}"
            else:
                raise ValueError(
                    "Input must be a list of fastq files or a single fastq file."
                )

        ref = self.snakemake.input.get("ref", "")
        if ref:
            self.extra = f" -r {ref}"

        fasta = self.snakemake.input.get("fasta", "")
        if fasta:
            self.extra = f" --genome-fasta-files {fasta}"

        fasta_dir = self.snakemake.input.get("fasta_dir", "")
        if fasta_dir:
            self.extra = f" --genome-fasta-directory {fasta_dir}"

        self.command = self.snakemake.params.get("command", "")

    def run(self):
        shell(
            "coverm"
            " {self.command}"
            " -t {self.snakemake.threads}"
            " {self.extra}"
            " -o {self.snakemake.output}"
            " {self.log}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)
