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


import tempfile
from os import path
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase
from snakemake_wrapper_utils.java import get_java_opts
from snakemake_wrapper_utils.samtools import get_samtools_opts


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        """Parser arguments
        """
        # Extract arguments.
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)
        sort = self.snakemake.params.get("sort", "none")
        self.sort_order = self.snakemake.params.get("sort_order", "coordinate")
        self.sort_extra = self.snakemake.params.get("sort_extra", "")
        self.samtools_opts = get_samtools_opts(
            self.snakemake, parse_threads=False, param_name="sort_extra"
        )
        self.java_opts = get_java_opts(self.snakemake)

        self.bwa_threads = self.snakemake.threads
        samtools_threads = self.snakemake.threads - 1

        index = self.snakemake.input.get("index", "")
        if isinstance(index, str):
            self.index = path.splitext(self.snakemake.input.idx)[0]
        else:
            self.index = path.splitext(self.snakemake.input.idx[0])[0]

        # Check inputs/arguments.
        if not isinstance(self.snakemake.input.reads, str) and len(self.snakemake.input.reads) not in {
            1,
            2,
        }:
            raise ValueError("input must have 1 (single-end) or 2 (paired-end) elements")

        if self.sort_order not in {"coordinate", "queryname"}:
            raise ValueError(f"Unexpected value for sort_order ({self.sort_order})")

        # Determine which pipe command to use for converting to bam or sorting.
        if sort == "none":
            # Correctly assign number of threads according to user request
            if samtools_threads >= 1:
                self.samtools_opts += f" --threads {samtools_threads} "

            if str(self.snakemake.output[0]).lower().endswith(("bam", "cram")):
                # Simply convert to bam using samtools view.
                self.pipe_cmd = " | samtools view {self.samtools_opts} > {self.snakemake.output[0]}"
            else:
                # Do not perform any sort nor compression, output raw sam
                self.pipe_cmd = " > {self.snakemake.output[0]} "

        elif sort == "samtools":
            # Correctly assign number of threads according to user request
            if samtools_threads >= 1:
                self.samtools_opts += f" --threads {samtools_threads} "

            # Add name flag if needed.
            if self.sort_order == "queryname":
                self.sort_extra += " -n"

            # Sort alignments using samtools sort.
            self.pipe_cmd = " | samtools sort {self.samtools_opts} {self.sort_extra} -T {tmpdir} > {self.snakemake.output[0]}"

        elif sort == "picard":
            # Correctly assign number of threads according to user request
            self.bwa_threads = self.bwa_threads - 1
            if self.bwa_threads <= 0:
                raise ValueError(
                    "Not enough threads requested. This wrapper requires exactly one more."
                )

            # Sort alignments using picard SortSam.
            self.pipe_cmd = (
                " | picard SortSam {self.java_opts} {self.sort_extra} "
                "--INPUT /dev/stdin "
                "--TMP_DIR {tmpdir} "
                "--SORT_ORDER {self.sort_order} "
                "--OUTPUT {self.snakemake.output[0]}"
            )

        else:
            raise ValueError(f"Unexpected value for params.sort ({sort})")

    def run(self):
        # run shell command
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "(bwa-mem2 mem"
                " -t {self.bwa_threads}"
                " {self.extra}"
                " {self.index}"
                " {self.snakemake.input.reads}"
                " " + self.pipe_cmd + " ) {self.log}"
            )
            shell("samtools index {self.snakemake.output[0]}")


if __name__ == "__main__":
    Wrapper(snakemake)