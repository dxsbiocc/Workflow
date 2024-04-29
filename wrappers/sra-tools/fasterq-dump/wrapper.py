# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/20 11:11:09
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import os
import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.snakemake import get_mem
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        # input
        input_acc = self.snakemake.input.get('sra', None)
        acc = self.snakemake.params.get("accession", None)
        assert input_acc ^ acc, "Please specify either input.sra or params.accession"
        self.accession = input_acc or acc
        # Parse memory
        mem_mb = get_mem(self.snakemake, "MiB")

        # Outdir
        outdir = os.path.dirname(self.snakemake.output[0])
        self.outdir = f"--outdir {outdir}" if outdir else ""

        # Output compression
        self.compress = ""
        mem = f"-m{mem_mb}" if mem_mb else ""

        for output in self.snakemake.output:
            out_name, out_ext = os.path.splitext(output)
            if out_ext == ".gz":
                self.compress += f"pigz -p {self.snakemake.threads} {out_name}; "
            elif out_ext == ".bz2":
                self.compress += f"pbzip2 -p{self.snakemake.threads} {mem} {out_name}; "

        self.mem = f"--mem {mem_mb}M" if mem_mb else ""


    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "(fasterq-dump"
                " {self.outdir}"
                " {self.extra}"
                " --temp {tmpdir} --threads {self.snakemake.threads} {self.mem}"
                " {self.accession};"
                " {self.compress}"
                ") {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)