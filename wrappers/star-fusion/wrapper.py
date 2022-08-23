# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/22 20:39:28
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
        # R1
        fq1 = self.snakemake.input.get("fastq1", None)
        self.left = f"--left_fq {fq1}" if fq1 else ""
        # R2
        fq2 = self.snakemake.input.get("fastq2", None)
        self.right = f"--left_fq {fq2}" if fq2 else ""
        # use star align result
        junc = self.snakemake.input.get("junction", None)
        self.junc = f"--chimeric_junction {junc}" if junc else ""
        if not (bool(fq1 or fq2) ^ bool(junc)):
            raise ValueError("You can only choose one from FASTQ file and STAR align junction file.")
        # index
        index = self.snakemake.input.get("index", None)
        assert index, "genome lib must be set."
        self.index = f"--genome_lib_dir {index}"


    def run(self):
        shell(
            "(STAR-Fusion"
            " {self.index}"
            " {self.left} {self.right}"
            " {self.junc}"
            " --output_dir {self.snakemake.output}"
            ") {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)