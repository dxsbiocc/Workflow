# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2024/06/22 14:11:34
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
        self.fasta = self.snakemake.input.get(
            "fasta"
        )  # Genomic FASTA file(s), comma-delimited
        self.gtf = self.snakemake.input.get(
            "gtf"
        )  # Reference GTF file(s), comma-delimited
        self.features = self.snakemake.input.get(
            "features"
        )  # [`kite` workflow only] Path to TSV containing barcodes and feature names

        self.mode = self.snakemake.params.get(
            "mode", "standard"
        )  # standard,nac,kite,custom
        if self.mode not in ["standard", "nac", "kite", "custom"]:
            raise ValueError("Invalid mode: {}".format(self.mode))

        # required arguments
        output = os.path.abspath(self.snakemake.output[0])
        if not os.path.exists(output):
            os.makedirs(output)
        self.index = os.path.join(
            output, "kb_ref.idx"
        )  # Path to the kallisto index to be constructed
        self.t2g = os.path.join(
            output, "t2g.txt"
        )  # Path to transcript-to-gene mapping to be generated
        self.f1 = os.path.join(
            output, "cdna.fa"
        )  # Path to the cDNA FASTA (lamanno, nucleus) or mismatch FASTA (kite) to be generated
        self.options = ""
        # required arguments for `lamanno` and `nucleus` workflows
        if self.mode == "nac":
            f2 = os.path.join(
                output, "intron.fa"
            )  # Path to the intron FASTA to be generated
            t2c = os.path.join(
                output, "cdna_t2c.txt"
            )  # Path to generate cDNA transcripts-to-capture
            t2i = os.path.join(
                output, "intron_t2c.txt"
            )  # Path to generate intron transcripts-to-capture
            self.options = f" -f2 {f2} -c1 {t2c} -c2 {t2i}"

    def run(self):
        shell(
            "kb ref "
            " -i {self.index}"
            " -g {self.t2g}"
            " -f1 {self.f1}"
            " {self.options}"
            " -t {self.snakemake.threads}"
            " {self.extra}"
            " --workflow {self.mode}"
            " {self.fasta}"  # positional arguments
            " {self.gtf}"
            " {self.features}"  # `kite` workflow only
            " {self.log}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)
