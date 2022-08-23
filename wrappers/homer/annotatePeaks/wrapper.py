# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/21 14:34:27
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import os
import sys
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.genome = self.snakemake.input.get("genome", "none")
        motif_files = self.snakemake.input.get("motif_files", None)
        matrix = self.snakemake.output.get("matrix", None)

        # optional files
        opt_files = {
            "gtf": "-gtf",
            "gene": "-gene",
            "motif_files": "-m",
            "filter_motiv": "-fm",
            "center": "-center",
            "nearest_peak": "-p",
            "tag": "-d",
            "vcf": "-vcf",
            "bed_graph": "-bedGraph",
            "wig": "-wig",
            "map": "-map",
            "cmp_genome": "-cmpGenome",
            "cmp_Liftover": "-cmpLiftover",
            "advanced_annotation": "-ann",
            "mfasta": "-mfasta",
            "mbed": "-mbed",
            "mlogic": "-mlogic",
        }

        requires_motives = False
        for opt in opt_files:
            file = None
            if opt in {"mfasta", "mbed", "mlogic"}:
                file = self.snakemake.output.get(opt, "")
                if file:
                    requires_motives = True
            else:
                file = self.snakemake.input.get(opt, "")
            if file:
                self.extra += " {flag} {file}".format(flag=opt_files[opt], file=file)

        if requires_motives and not motif_files:
            sys.exit(
                "The optional output files require motif_file(s) as input. For more information please see http://homer.ucsd.edu/homer/ngs/annotation.html."
            )

        # optional matrix output files:
        if matrix:
            if not motif_files:
                sys.exit(
                    "The matrix output files require motif_file(s) as input. For more information please see http://homer.ucsd.edu/homer/ngs/annotation.html."
                )
            ext = ".count.matrix.txt"
            matrix_out = [i for i in self.snakemake.output if i.endswith(ext)][0]
            matrix_name = os.path.basename(matrix_out[: -len(ext)])
            self.extra += " -matrix {}".format(matrix_name)

    def run(self):
        shell(
            "(annotatePeaks.pl"
            " {self.snakemake.params.mode}"
            " {self.snakemake.input.peaks}"
            " {self.genome}"
            " {self.extra}"
            " -cpu {self.snakemake.threads}"
            " > {self.snakemake.output.annotations}"
            ") {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)