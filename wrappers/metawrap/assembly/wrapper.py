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
        r1 = self.snakemake.input.get("r1", "")
        if isinstance(r1, str):
            self.r1 = [r1]
        r2 = self.snakemake.input.get("r2", "")
        if isinstance(r2, str):
            self.r2 = [r2]

        self.output = os.path.basename(str(self.snakemake.output))

        assembly_method = self.snakemake.params.get("assembly_method", "megahit")
        if assembly_method == "megahit":
            self.extra += " --megahit"
        elif assembly_method == "metaspades":
            self.extra += " --metaspades"
        else:
            raise ValueError(f"Invalid assembly method: {assembly_method}")

        self.memory_requirements = get_mem(self.snakemake, out_unit="MiB")

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            # ln -s the input files to the temporary directory
            new_r1s = []
            new_r2s = []
            for r1 in self.r1:
                new_r1 = r1.replace(".R1.fq.gz", "_1.fastq")
                os.symlink(r1, os.path.join(tmpdir, os.path.basename(new_r1)))
                new_r1s.append(new_r1)
            for r2 in self.r2:
                new_r2 = r2.replace(".R2.fq.gz", "_2.fastq")
                os.symlink(r2, os.path.join(tmpdir, os.path.basename(new_r2)))
                new_r2s.append(new_r2)
            shell(
                "metawrap assembly"
                " -1 {new_r1s} -2 {new_r2s}"
                " -o {self.output}"
                " -t {self.snakemake.threads}"
                " -m {self.memory_requirements}"
                " {self.extra}"
                " {self.log}"
            )


if __name__ == "__main__":
    Wrapper(snakemake)
