# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/28 10:14:30
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import tempfile
from pathlib import Path
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        pass
            

    def run(self):
        if Path(self.snakemake.input.bed).stat().st_size:
            with tempfile.TemporaryDirectory() as tmpdir:
                shell(
                    "get_seqs {self.extra}"
                    " -p {tmpdir}/out"
                    " {self.snakemake.input.bed}"
                    " {self.snakemake.input.fas}"
                    " {self.log}"
                )

                if self.snakemake.output.get("hap"):
                    shell("cat {tmpdir}/out.hap.fa > {self.snakemake.output.hap}")
                if self.snakemake.output.get("purged"):
                    shell("cat {tmpdir}/out.purged.fa > {self.snakemake.output.purged}")
        else:
            # If BED file empty, copy input to output since `get_seqs` will segfault
            log = Path(self.snakemake.log[0])
            log.write_text(
                "WARN: Input BED file is empty. Input FASTA file will be copied to output."
            )
            shell("cp {self.snakemake.input.fas} {self.snakemake.output.hap}")
            Path(self.snakemake.output.purged).touch()


if __name__ == '__main__':
    Wrapper(snakemake)