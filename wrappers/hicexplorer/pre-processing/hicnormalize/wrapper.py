# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/31 20:41:02
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
        self.normalize = self.snakemake.params.get("normalize", "smallest")
        assert self.normalize in [
            "smallest", 
            "norm_range", 
            "multiplicative"
        ], "normalize must be one of [smallest, norm_range, multiplicative]"

        assert len(self.snakemake.input) == len(self.snakemake.output), "input and output must be the same length"

        if self.normalize == 'multiplicative':
            self.extra += f" --multiplicativeValue {self.snakemake.params.get('multiplicativeValue', 1)}"

        self.setToZeroThreshold = self.snakemake.params.get("setToZeroThreshold", 0.0)
        if self.normalize == 'smallest' and self.setToZeroThreshold != 1:
            Warning("It is recommended to set it for the normalize mode “smallest” to 1.0!")
    
    def run(self):
        shell(
            "(hicNormalize"
            " {self.extra}"
            " --setToZeroThreshold {self.setToZeroThreshold}"
            " --matrices {self.snakemake.input}"
            " --normalize {self.normalize}"
            " --outFileName {self.snakemake.output}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)