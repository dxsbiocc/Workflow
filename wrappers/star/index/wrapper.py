# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/20 21:17:54
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import tempfile
from snakemake.shell import shell
from snakemake.utils import makedirs
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.sjdb_overhang = self.snakemake.params.get("sjdbOverhang", "100")

        gtf = self.snakemake.input.get("gtf")
        if gtf is not None:
            self.gtf = f"--sjdbGTFfile {gtf}"
            self.sjdb_overhang = f"--sjdbOverhang {self.sjdb_overhang}"
        else:
            self.gtf, self.sjdb_overhang = "", ""

        makedirs(self.snakemake.output)

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "STAR"
                " --runThreadN {self.snakemake.threads}"  # Number of threads
                " --runMode genomeGenerate"  # Indexation mode
                " --genomeDir {self.snakemake.output}"  # Path to output
                " --genomeFastaFiles {self.snakemake.input.fasta}"  # Path to fasta files
                " {self.sjdb_overhang}"  # Read-len - 1
                " {self.gtf}"  # Highly recommended GTF
                " {self.extra}"  # Optional parameters
                " --outTmpDir {tmpdir}/STARtmp"  # Temp dir
                " {self.log}"  # Logging
            )


if __name__ == '__main__':
    Wrapper(snakemake)