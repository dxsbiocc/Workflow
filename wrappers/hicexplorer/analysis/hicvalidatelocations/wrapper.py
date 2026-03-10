# -*- encoding: utf-8 -*-
# ============================================================
# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/09/04 15:09:36
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import os
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):
    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.validtype = self.snakemake.params.get("validtype", "bed")
        assert self.validtype in ["bed", "cool"], "validtype must be in ['bed', 'cool']"

        self.method = self.snakemake.params.get("method", "loops")
        assert self.method in ["loops", "tad"], "method must be in ['loops', 'tad']"

        self.outFileName = os.path.commonprefix(self.snakemake.output).strip("_")
    
    def run(self):
        shell(
            "(hicValidateLocations"
            " {self.extra}"
            " --data {self.snakemake.input.data}"
            " --validationData {self.snakemake.input.valid}"
            " --validationType {self.validtype}"
            " --method {self.method}"
            " --resolution {self.snakemake.params.resolution}"
            " --outFileName {self.outFileName}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)