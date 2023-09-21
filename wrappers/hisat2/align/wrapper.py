# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/08 21:35:57
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
        # Input file wrangling
        reads = self.snakemake.input.get("reads")
        if isinstance(reads, str):
            self.input_flags = "-U {0}".format(reads)
        elif len(reads) == 1:
            self.input_flags = "-U {0}".format(reads[0])
        elif len(reads) == 2:
            self.input_flags = "-1 {0} -2 {1}".format(*reads)
        else:
            raise RuntimeError(
                "Reads parameter must contain at least 1 and at most 2" " input files."
            )
        # index prefix
        self.index = os.path.commonprefix(self.snakemake.input.get("index")).strip('.')
        # sort bam
        sort = self.snakemake.params.get("sort", "samtools")
        sort_order = self.snakemake.params.get("sort_order", "coordinate")
        sort_extra = self.snakemake.params.get("sort_extra", "")

        if sort_order not in {"coordinate", "queryname"}:
            raise ValueError("Unexpected value for sort_order ({})".format(sort_order))

        # Determine which pipe command to use for converting to bam or sorting.
        if sort == "none":
            # Simply convert to bam using samtools view.
            self.pipe_cmd = f"samtools view -Sbh -o {self.snakemake.output[0]} -"
        elif sort == "samtools":
            # Add name flag if needed.
            if sort_order == "queryname":
                sort_extra += " -n"
            # Sort alignments using samtools sort.
            self.pipe_cmd = f"samtools sort -T {self.tmp} {sort_extra} -o {self.snakemake.output[0]} -"
        elif sort == "picard":
            # Sort alignments using picard SortSam.
            self.pipe_cmd = (
                f"picard SortSam {sort_extra} --INPUT /dev/stdin"
                f" --OUTPUT {self.snakemake.output[0]} --SORT_ORDER {sort_order} --TMP_DIR {self.tmp}"
            )
        else:
            raise ValueError("Unexpected value for params.sort ({})".format(sort))

    def run(self):
        # Executed shell command
        shell(
            "(hisat2 {self.extra}"
            " --threads {self.snakemake.threads}"
            " -x {self.index}"
            " {self.input_flags}"
            " | {self.pipe_cmd})"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)