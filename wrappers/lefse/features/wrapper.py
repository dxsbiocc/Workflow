# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 09:29:47
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file
# ============================================================


import os
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.dataset = self.snakemake.input.dataset
        self.lefse = self.snakemake.input.lefse
        self.output_dir = str(self.snakemake.output)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    def run(self):
        shell(
            "lefse_plot_features.py {self.extra} {self.dataset} {self.lefse} {self.output_dir}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)
