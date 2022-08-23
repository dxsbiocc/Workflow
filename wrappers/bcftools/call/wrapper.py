# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/21 15:25:37
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
        self.bcftools_opts = get_bcftools_opts(self.snakemake, parse_ref=False, parse_memory=False)

        valid_caller_opts = {"-c", "--consensus-caller", "-m", "--multiallelic-caller"}

        self.caller_opt = self.snakemake.params.get("caller", "")
        if self.caller_opt.strip() not in valid_caller_opts:
            raise ValueError(
                "bcftools call expects either -m/--multiallelic-caller or "
                "-c/--consensus-caller as caller option."
            )

    def run(self):
        shell(
            "bcftools call"
            " {self.bcftools_opts}"
            " {self.caller_opt}"
            " {self.extra}"
            " {self.snakemake.input[0]}"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)