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
        reads = self.snakemake.input.get("reads")
        # use star align result
        junc = self.snakemake.input.get("junction")
        if reads:
            if len(reads) == 2:
                fq1, fq2 = reads
                self.input = f"--left_fq {fq1} --right_fq {fq2}"
            else:
                fq1 = reads[0] if isinstance(reads, list) else reads
                self.input = f"--left_fq {fq1}"
        elif junc:
            self.input = f"--chimeric_junction {junc}" if junc else ""
        else:
            raise ValueError("Neither 'junction' or 'reads' not exist")
        
        if not (bool(reads) ^ bool(junc)):
            raise ValueError("You can only choose one from FASTQ file and STAR align junction file.")
        # index
        index = self.snakemake.input.get("index", None)
        assert index, "genome lib must be set."
        self.index = f"--genome_lib_dir {index}"


    def run(self):
        shell(
            "(STAR-Fusion"
            " {self.index}"
            " {self.input}"
            " --output_dir {self.snakemake.output}"
            ") {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)