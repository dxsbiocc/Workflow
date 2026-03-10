# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/25 21:37:48
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import os
import glob
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):
    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.wd = os.path.abspath(self.snakemake.input[0])

        run_type = "[generic|specific]*"
        file_list = glob.glob(f'{self.wd}/short_summary.{run_type}.*.*.txt')
        if len(file_list) == 0:
            raise Exception("No short_summary file found! please check your file format.\
                            The file name should be short_summary.[generic|specific].[dataset].[edit_name_here].txt")

    def run(self):
        shell(
            "generate_plot.py"
            " -wd {self.wd}"
            " {self.extra}"
            " {self.log}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)