# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/31 17:42:36
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
        ## Extract arguments
        expected = self.snakemake.input.get("expected", "")
        if expected:
            self.extra += f"--expected {expected}"

        view = self.snakemake.input.get("view", "")
        if view:
            self.extra = f"--view {view}"

        self.resolution = self.snakemake.params.get(
            "resolution", self.snakemake.wildcards.get("resolution", 0)
        )
        if not self.resolution:
            raise ValueError("Please specify resolution either as a wildcard or as a parameter")


    def run(self):
        shell(
            "(cooltools pileup"
            " {self.extra}"
            " {self.snakemake.input.cooler}::resolutions/{self.resolution}"
            " {self.snakemake.input.features}"
            " --features-format {self.snakemake.params.features_format}"
            " -p {self.snakemake.threads}"
            " -o {self.snakemake.output}"
            ") {self.log}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)