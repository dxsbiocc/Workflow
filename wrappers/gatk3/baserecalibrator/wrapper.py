# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/09 20:48:15
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

        self.known = self.snakemake.input.get("known", "")
        if self.known:
            if isinstance(self.known, str):
                self.known = [self.known]
            self.known = list(map("--knownSites {}".format, self.known))

        bed = self.snakemake.params.get("bed", "")
        self.bed = f"--intervals {bed}" if bed else ""
            
    def run(self):
        shell(
            "gatk3 {self.java_opts}"
            " --analysis_type BaseRecalibrator"
            " --num_cpu_threads_per_data_thread {self.snakemake.threads}"
            " --input_file {self.snakemake.input.bam}"
            " {self.known}"
            " --reference_sequence {self.snakemake.input.ref}"
            " {self.bed}"
            " {self.extra}"
            " --out {self.snakemake.output}"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)