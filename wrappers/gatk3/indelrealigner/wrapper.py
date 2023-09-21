# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/09 20:47:23
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.java_opts = get_java_opts(self.snakemake)

        bed = self.snakemake.params.get("bed", "")
        self.bed = f"--intervals {bed}" if bed else ""
        
        self.known = self.snakemake.input.get("known", "")
        if self.known:
            if isinstance(self.known, str):
                self.known = f"--knownAlleles {self.known}"
            else:
                self.known = list(map("----knownAlleles {}".format, self.known))

        output_bai = self.snakemake.output.get("bai", None)
        if output_bai is None:
            self.extra += " --disable_bam_indexing"

    def run(self):
        shell(
            "gatk3 {self.java_opts}"
            " --analysis_type IndelRealigner"
            " --input_file {self.snakemake.input.bam}"
            " --reference_sequence {self.snakemake.input.ref}"
            " {self.known}"
            " {self.bed}"
            " --targetIntervals {self.snakemake.input.target_intervals}"
            " {self.extra}"
            " --out {self.snakemake.output.bam}"
            " {self.log}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)