# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/18 15:29:35
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
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)

        n = len(self.snakemake.input)
        self.out_file = None
        # data type
        if n == 1:   # single-end
            self.out_file = f" -o {self.snakemake.output.trimmed[0]}"
        elif n == 2:  # paired-end
            self.out_file = f" -o {self.snakemake.output.trimmed[0]} -p {self.snakemake.output.trimmed[1]}"
        else:
            raise ValueError("Input reads must contain 1 (single-end) or 2 (paired-end) elements.")

        self.adapters = self.snakemake.params.get("adapters", "")
        
        assert (
            self.extra != "" or self.adapters != ""
        ), "No options provided to cutadapt. Please use 'params: adapters=' or 'params: extra='."

    def run(self):
        shell(
            "cutadapt"
            " --cores {self.snakemake.threads}"
            " {self.adapters}"
            " {self.extra}"
            " {self.out_file}"
            " {self.snakemake.input}"
            " > {self.snakemake.output.qc} {self.log}"
        )

if __name__ == '__main__':
    Wrapper(snakemake)