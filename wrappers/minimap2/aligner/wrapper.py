# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/08 21:17:09
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from os import path
from snakemake.shell import shell
from snakemake_wrapper_utils.samtools import infer_out_format
from snakemake_wrapper_utils.samtools import get_samtools_opts
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        samtools_opts = get_samtools_opts(self.snakemake, parse_output=False)
        sort = self.snakemake.params.get("sorting", "none")
        sort_extra = self.snakemake.params.get("sort_extra", "")
        out_ext = infer_out_format(self.snakemake.output[0])

        self.pipe_cmd = ""
        if out_ext != "PAF":
            # Add option for SAM output
            self.extra += " -a"
            # Determine which pipe command to use for converting to bam or sorting.
            if sort == "none":
                if out_ext != "SAM":
                    # Simply convert to output format using samtools view.
                    self.pipe_cmd = f"| samtools view -h {samtools_opts}"
            elif sort in ["coordinate", "queryname"]:
                # Add name flag if needed.
                if sort == "queryname":
                    sort_extra += " -n"
                # Sort alignments.
                self.pipe_cmd = f"| samtools sort {sort_extra} {samtools_opts}"
            else:
                raise ValueError(f"Unexpected value for params.sort: {sort}")
        if out_ext == "GZ":
            self.pipe_cmd += " | pigz -c - "

    def run(self):
        shell(
            "(minimap2"
            " -t {self.snakemake.threads}"
            " {self.extra} "
            " {self.snakemake.input.target}"
            " {self.snakemake.input.query}"
            " {self.pipe_cmd}"
            " > {self.snakemake.output[0]}"
            ") {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)