# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/04/21 14:49:21
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import os

from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        # the reference_name argument is inferred by stripping the .seq suffix from
        # the output.seq value
        output_directory = os.path.dirname(os.path.abspath(self.snakemake.output.seq))
        seq_file = os.path.basename(self.snakemake.output.seq)
        if seq_file.endswith(".seq"):
            self.reference_name = os.path.join(output_directory, seq_file[:-4])
        else:
            raise Exception("output.seq has an invalid file suffix (must be .seq)")

        for output_variable, output_path in self.snakemake.output.items():
            if not os.path.abspath(output_path).startswith(self.reference_name):
                raise Exception(
                    "the path for {} is inconsistent with that of output.seq".format(output_variable)
                )
            
        gtf = self.snakemake.input.get("gtf", None)
        if gtf:
            self.extra += f" --gtf {gtf}"

        mapping = self.snakemake.params.get("mapping")
        if mapping == "bowtie2":
            self.extra += " --bowtie2"
        elif mapping == "star":
            self.extra += " --star"
        elif mapping == "hisat2":
            self.extra += " --hisat2-hca"

    def run(self):
        shell(
            "rsem-prepare-reference"
            " --num-threads {self.snakemake.threads}"
            " {self.extra}"
            " {self.snakemake.input.reference_genome}"
            " {self.reference_name} "
            "{self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)