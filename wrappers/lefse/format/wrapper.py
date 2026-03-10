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
        self.input = self.snakemake.input
        self.lefse_input = self.snakemake.output.lefse_input
        self.lefse_format = self.snakemake.output.lefse_format

        self.group_line = self.snakemake.params.group_line

    def run(self):
        # generate lefse input file
        with open(str(self.input), "r") as f:
            lines = f.readlines()
        lines[0] = lines[0].replace("ID", "subject_id", 1)
        site_line = lines[0].replace("subject_id", "site", 1)
        lines.insert(0, site_line)
        lines.insert(0, self.group_line)
        # create lefse directory
        if not os.path.exists(os.path.dirname(self.lefse_input)):
            os.makedirs(os.path.dirname(self.lefse_input))
        # write lefse input file
        with open(self.lefse_input, "w") as f:
            f.writelines(lines)

        if not os.path.exists(self.lefse_input):
            raise FileNotFoundError(f"Lefse input file {self.lefse_input} not found")

        shell(
            "lefse_format_input.py {self.lefse_input} {self.lefse_format} {self.extra}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)
