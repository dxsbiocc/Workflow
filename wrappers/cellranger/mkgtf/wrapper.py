# -*- encoding: utf-8 -*-
# ============================================================
# File        : cellranger_count.py
# Time        : 2024/06/20 11:17:39
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

    def __init__(self, snakemake):
        super().__init__(snakemake)

    def parser(self):

        self.input_gtf = self.snakemake.input.gtf
        self.output_gtf = self.snakemake.output
        
        self.cellranger_path = self.snakemake.params.get('cellranger')
        if not os.path.exists(self.cellranger_path):
            raise ValueError("cellranger path not exists")
        self.genome = self.snakemake.params.get('genome')
    
    def run(self):
        shell(
            "export PATH=PATH:{self.cellranger_path}/bin &&"
            " {self.cellranger_path}/bin/cellranger mkgtf "
            " {self.input_gtf} "
            " {self.output_gtf} "
            " {self.extra}"
            " {self.log}"
        )

if __name__ == '__main__':
    wrapper = Wrapper(snakemake)