# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/20 11:45:54
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
        # list of accessions
        self.accessions = self.snakemake.params.get("accessions", "")
        # accession file
        infile = self.snakemake.input.get("infile", None)
        # checkout input
        if not (bool(self.accessions) ^ bool(infile)):
            raise ValueError("accessions or infile must choose one!")
        self.infile = f" --option-file {infile}" if infile else "" if infile else ""

        # checkout output
        self.output_file, self.output_dir = "", ""
        out_file = self.snakemake.output.get('output_file', None)
        if out_file:
            if isinstance(self.accessions, list) or infile:
                raise ValueError("Download multiple accessions, but set a single file!")
            self.output_file = f"--output-file {out_file}"
        out_dir = self.snakemake.output.get('output_dir', None)
        if out_dir:
            self.output_dir = f"--output-directory {out_dir}"


    def run(self):
        shell("(prefetch"
        " {self.accessions}"
        " {self.infile}"
        " {self.extra}"
        " {self.output_file}"
        " {self.output_dir}"
        ") {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)