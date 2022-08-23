# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/21 16:00:59
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase
from snakemake_wrapper_utils.bcftools import get_bcftools_opts


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.bcftools_opts = get_bcftools_opts(
            self.snakemake, parse_ref=False, parse_output_format=False, parse_memory=False
        )

        if "--tbi" in self.extra or "--csi" in self.extra:
            raise ValueError(
                "You have specified index format (`--tbi/--csi`) in `params.extra`; this is automatically infered from the first output file."
            )

        if self.snakemake.output[0].endswith(".tbi"):
            self.extra += " --tbi"
        elif self.snakemake.output[0].endswith(".csi"):
            self.extra += " --csi"
        else:
            raise ValueError("invalid index file format ('.tbi', '.csi').")

    def run(self):
        shell(
            "(bcftools index"
            " {self.bcftools_opts}"
            " {self.extra}"
            " {self.snakemake.input[0]}"
            ") {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)