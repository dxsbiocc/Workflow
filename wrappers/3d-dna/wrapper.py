# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/09/04 21:54:03
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
        self.out_dir = os.path.abspath(self.snakemake.output[0])
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

        self.mode = self.snakemake.params.get('mode', 'haploid')
        assert self.mode in ['haploid', 'diploid'], "mode must be one of ['haploid', 'diploid']"

        input_size = self.snakemake.params.get('input_size', 15000)
        if input_size:
            self.extra += f" -i {input_size}"
        rounds = self.snakemake.params.get('rounds', 2)
        if rounds:
            self.extra += f" -r {rounds}"
        stage = self.snakemake.params.get('stage', '')
        if stage:
            self.extra += f" -s {stage}"

        self.fasta = os.path.abspath(self.snakemake.input.get('fasta'))
        self.mnd = os.path.abspath(self.snakemake.input.get('mnd'))

        self.log = f'> {os.path.abspath(self.snakemake.log[0])} 2>&1'
    
    def run(self):
        shell(
            " cd {self.out_dir} &&"
            " (3d-dna "
            " --mode {self.mode}"
            " {self.fasta}"
            " {self.mnd}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)