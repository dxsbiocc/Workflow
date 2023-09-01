# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/31 17:33:59
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)
        ## Extract arguments
        window = self.snakemake.params.get("window", "")
        if isinstance(window, list):
            self.window = " ".join([str(w) for w in window])
        else:
            self.window = str(window)

        view = self.snakemake.input.get("view", "")
        if view:
            self.extra += f"--view {view}"
        self.chunksize = self.snakemake.params.get("chunksize", 20000000)

        self.resolution = self.snakemake.params.get(
            "resolution", self.snakemake.wildcards.get("resolution", 0)
        )
        if not self.resolution:
            raise ValueError("Please specify resolution either as a wildcard or as a parameter")


    def run(self):
        shell(
            "(cooltools insulation"
            " {self.extra} "
            " {self.snakemake.input.cooler}::resolutions/{self.resolution} "
            " {self.window}"
            " --chunksize {self.chunksize} "
            " -p {self.snakemake.threads} "
            " -o {self.snakemake.output}"
            ") {self.log}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)