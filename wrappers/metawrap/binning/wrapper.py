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


import os
import tempfile
from snakemake import shell
from snakemake_wrapper_utils.base import WrapperBase
from snakemake_wrapper_utils.snakemake import get_mem


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.assembly = self.snakemake.input.get("assembly", "")
        self.fastqs = self.snakemake.input.get("fastq", "")
        self.output = self.snakemake.output

        self.memory_requirements = get_mem(self.snakemake, out_unit="MiB")

        methods = self.snakemake.params.get("methods", "")
        if isinstance(methods, list):
            self.methods = " --".join(methods[:3])  # only take the first 3 methods
        else:
            self.methods = methods

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            new_fastqs = []
            for fastq in self.fastqs:
                new_fastq = (
                    os.path.basename(fastq)
                    .replace(".R1.fq.gz", "_1.fastq")
                    .replace(".R2.fq.gz", "_2.fastq")
                )
                shell(f"ln -s {fastq} {os.path.join(tmpdir, new_fastq)}")
                new_fastqs.append(os.path.join(tmpdir, new_fastq))
            shell(
                "metawrap binning"
                " -a {self.assembly}"
                " -o {self.output}"
                " -t {self.snakemake.threads}"
                " -m {self.memory_requirements}"
                " --{self.methods}"
                " {self.extra}"
                " {new_fastqs}"
                " {self.log}"
            )


if __name__ == "__main__":
    Wrapper(snakemake)
