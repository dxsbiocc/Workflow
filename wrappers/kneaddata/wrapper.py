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


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        # input files
        fq_list = self.snakemake.input.get("fq_list")
        if isinstance(fq_list, list) and len(fq_list) == 2:
            self.fastq = f"--input1 {fq_list[0]} --input2 {fq_list[1]}"
        elif isinstance(fq_list, list) and len(fq_list) == 1:
            self.fastq = f"--unpaired {fq_list[0]}"
        else:
            raise ValueError("input->fq_list must be a list of 1 or 2 elements.")

        # host genome database
        self.db = f"--reference-db {self.snakemake.input.get("db")}"
        # output files
        self.outdir = os.path.dirname(self.snakemake.output)
        # extra parameters
        self.extra = self.snakemake.params.get("extra")

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "kneaddata"
                " {self.fastq}"
                " {self.db}"
                " {self.extra}"
                " -t {self.snakemake.threads}"
                " --output {tmpdir}"
                " {self.log}"
            )

            # Summary of Quality Control Results (质控结果汇总)
            shell(
                "kneaddata_read_count_table"
                " --input {tmpdir}"
                " --output {self.outdir}/kneaddata.txt"
            )

            if os.path.exists(os.path.join(tmpdir, "*_kneaddata.fastq")):
                shell("mv {tmpdir}/*_kneaddata.fastq {self.snakemake.output.fastq}")

            if os.path.exists(os.path.join(tmpdir, "*_kneaddata.log")):
                shell("mv {tmpdir}/*_kneaddata.log {self.log}")

            if os.path.exists(os.path.join(tmpdir, "*contam.fastq")):
                shell("mv {tmpdir}/*contam.fastq {self.snakemake.output.contam}")

        os.remove(tmpdir)


if __name__ == "__main__":
    Wrapper(snakemake)
