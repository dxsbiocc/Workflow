# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2024/01/11 20:01:59
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file
# ============================================================


from tempfile import TemporaryDirectory
from snakemake import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.input = self.snakemake.input.fasta
        self.output = self.snakemake.output.cluster

    def run(self):
        with TemporaryDirectory() as tmpdir:
            shell(
                "mmseqs easy-cluster"
                " {self.input}"
                " {tmpdir}/easy-cluster_result"
                " {tmpdir}"
                " {self.extra}"
                " --threads {self.snakemake.threads}"
                " {self.log}"
            )
            # mv
            shell("mv {tmpdir}/easy-cluster_result_rep_seq.fasta {self.output}")


if __name__ == "__main__":
    Wrapper(snakemake)
