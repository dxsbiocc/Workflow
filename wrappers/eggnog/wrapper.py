# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/18 20:07:45
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file
# ============================================================


import os
import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.data_dir = self.snakemake.input.get("data_dir", "")
        self.input = self.snakemake.input.get("protein", "")
        self.output_dir = os.path.dirname(str(self.snakemake.output))
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        self.prefix = self.snakemake.params.get("prefix", "")

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "emapper.py"
                " --data_dir {self.data_dir}"
                " -i {self.input} "
                " --output_dir {self.output_dir} "
                " --output {self.prefix} "
                " {self.extra} "
                " {self.log}"
            )


if __name__ == "__main__":
    Wrapper(snakemake)
