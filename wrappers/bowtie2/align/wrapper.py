# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/18 16:33:36
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import os
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        # reads
        n = len(self.snakemake.input.reads)
        if n == 1:
            self.reads = "-U {}".format(*self.snakemake.input.reads)
        elif n == 2:
            self.reads = "-1 {} -2 {}".format(*self.snakemake.input.reads)
        else:
            raise ValueError("input->sample must have 1 (single-end) or 2 (paired-end) elements.")

        # index
        index = self.snakemake.input.get("index", None)
        assert index, "`index` parameter not set in `input`."
        if isinstance(index, str):
            self.index = os.path.splitext(index)[0]
        else:
            # the largest common prefix
            self.index = os.path.commonprefix(index).rstrip(".")
        # Extract arguments.
        sort = self.snakemake.params.get("sort", "none")
        sort_order = self.snakemake.params.get("sort_order", "coordinate")
        sort_extra = self.snakemake.params.get("sort_extra", "")
        # Determine which pipe command to use for converting to bam or sorting.
        if sort == "none":
            # Simply convert to bam using samtools view.
            self.pipe_cmd = f"samtools view -Sbh -o {self.snakemake.output[0]}"
        elif sort == "samtools":
            # Add name flag if needed.
            if sort_order == "queryname":
                sort_extra += " -n"
            # Sort alignments using samtools sort.
            self.pipe_cmd = f"samtools sort {sort_extra} -o {self.snakemake.output[0]}"
        elif sort == "picard":
            # Sort alignments using picard SortSam.
            self.pipe_cmd = (
                f"picard SortSam {sort_extra} --INPUT /dev/stdin"
                f" --OUTPUT {self.snakemake.output[0]} --SORT_ORDER {sort_order} --TMP_DIR {self.tmp}"
            )
        else:
            raise ValueError("Unexpected value for params.sort ({})".format(sort))

    def run(self):
        shell(
            "(bowtie2"
            " --threads {self.snakemake.threads}"
            " {self.reads} "
            " -x {self.index}"
            " {self.extra} "
            "| " + self.pipe_cmd + ") {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)