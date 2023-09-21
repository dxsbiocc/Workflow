# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/04/21 14:12:19
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
        bam = self.snakemake.input.get("bam", "")
        fq_one = self.snakemake.input.get("fq_one", "")
        fq_two = self.snakemake.input.get("fq_two", "")

        if bam:
            if fq_one:
                raise Exception("Only input.bam or input.fq_one expected, got both.")
            self.input_bam = "--alignments"
            self.input_string = bam
            paired_end = self.snakemake.params.get("paired_end", False)
        else:
            self.input_bam = ""
            if fq_one:
                if isinstance(fq_one, list):
                    num_fq_one = len(fq_one)
                    input_string = ",".join(fq_one)
                else:
                    num_fq_one = 1
                    input_string = fq_one
                if fq_two:
                    paired_end = True
                    if isinstance(fq_two, list):
                        num_fq_two = len(fq_two)
                        if num_fq_one != num_fq_two:
                            raise Exception(
                                "Got {} R1 FASTQs, {} R2 FASTQs.".format(num_fq_one, num_fq_two)
                            )
                    else:
                        fq_two = [fq_two]
                    input_string += " " + ",".join(fq_two)
                else:
                    paired_end = False
            else:
                raise Exception("Expected input.bam or input.fq_one, got neither.")
        if paired_end:
            self.paired_end_string = "--paired-end"
        else:
            self.paired_end_string = ""

        genes_results = self.snakemake.output.genes_results
        if genes_results.endswith(".genes.results"):
            self.output_prefix = genes_results[: -len(".genes.results")]
        else:
            raise Exception(
                "output.genes_results file name malformed "
                "(rsem will append .genes.results suffix)"
            )
        if not self.snakemake.output.isoforms_results.endswith(".isoforms.results"):
            raise Exception(
                "output.isoforms_results file name malformed "
                "(rsem will append .isoforms.results suffix)"
            )

        self.reference_prefix = os.path.splitext(self.snakemake.input.reference[0])[0]

    def run(self):
        shell(
            "rsem-calculate-expression"
            " --num-threads {self.snakemake.threads}"
            " {self.extra}"
            " {self.paired_end_string}"
            " {self.input_bam}"
            " {self.input_string}"
            " {self.reference_prefix}"
            " {self.output_prefix}"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)