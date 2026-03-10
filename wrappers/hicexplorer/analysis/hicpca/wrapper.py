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


from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):
    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        output_fmt = self.snakemake.params.get("output_fmt", "")
        assert output_fmt in ["bedgraph", "bigwig"], "'output_fmt' must be one of bedgraph, bigwig"
        if output_fmt:
            self.extra += f" --format {output_fmt}"

        method = self.snakemake.params.get("method", "")
        assert method in ["dist_norm", "lieberman"], "'method' must be one of dist_norm, lieberman"
        if method:
            self.extra += f" --method {method}"

        if we := self.snakemake.params.get("we", ""):
            self.extra += f" --whichEigenvectors {we}"

    def run(self):
        shell(
            "(hicPCA"
            " {self.extra}"
            " --matrix {self.snakemake.input}"
            " --outputFileName {self.snakemake.output}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)