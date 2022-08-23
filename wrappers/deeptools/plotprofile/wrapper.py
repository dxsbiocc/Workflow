# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 15:53:52
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
        out_region = self.snakemake.output.get("regions")
        out_data = self.snakemake.output.get("data")

        self.optional_output = ""
        if out_region:
            self.optional_output += f" --outFileSortedRegions {out_region} "

        if out_data:
            self.optional_output += f" --outFileNameData {out_data} "

    def run(self):
        shell(
            "(plotProfile "
            "-m {self.snakemake.input[0]} "
            "-o {self.snakemake.output.plot_img} "
            "{self.optional_output} "
            "{self.extra}) {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)