# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/17 15:55:53
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import os
import re
import sys
import tempfile

from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)
        # Extract arguments.
        sort = self.snakemake.params.get("sort", "samtools")
        sort_order = self.snakemake.params.get("sort_order", "coordinate")
        sort_extra = self.snakemake.params.get("sort_extra", "")

        # get reference genome file
        index = self.snakemake.input.get("index", None)
        assert index, "`index` parameter not set in `input`."

        if isinstance(index, str):
            self.index = os.path.splitext(index)[0]
        else:  
            # the largest common prefix
            self.index = os.path.commonprefix(index).rstrip(".")

        # tmp dir
        if re.search(r"-T\b", sort_extra) or re.search(r"--TMP_DIR\b", sort_extra):
            sys.exit(
                "You have specified temp dir (`-T` or `--TMP_DIR`) in params.sort_extra; this is automatically set from params.tmp_dir."
            )

        # Check inputs/arguments.
        if not isinstance(self.snakemake.input.reads, str) and len(self.snakemake.input.reads) not in {
            1,
            2,
        }:
            raise ValueError("input must have 1 (single-end) or " "2 (paired-end) elements")

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
        # run shell command
        with tempfile.TemporaryDirectory() as tmp:
            shell(
                "(bwa mem"
                " -t {self.snakemake.threads}"
                " {self.extra}"
                " {self.index}"
                " {self.snakemake.input.reads}"
                " | {self.pipe_cmd})"
                " {self.log}"
            )
            shell("samtools index {self.snakemake.output[0]}")


if __name__ == "__main__":
    Wrapper(snakemake)
