# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/25 16:48:53
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
        self.db = os.path.abspath(self.snakemake.input.get("db"))
        hifi = self.snakemake.input.get("hifi")
        if isinstance(hifi, list):
            self.hifi = " ".join([os.path.abspath(i) for i in hifi])
        else:
            self.hifi = os.path.abspath(hifi)

        self.output = self.snakemake.output[0]
        if not os.path.exists(self.output):
            os.makedirs(self.output)

        self.output_prefix = self.snakemake.params.get("output_prefix", "merqury")

    def run(self):
        shell(
            "(cd {self.output} &&"
            " merqury.sh "
            " {self.db}"
            " {self.hifi}"
            " {self.output_prefix}"
            " {self.extra}"
            " ) {self.log}"
        )
        # remove 
        # log_dir = os.path.dirname(self.log).strip('> ')
        # shell('mv {self.output}/logs/* {log_dir}')
        # shell('rm -rf {self.output}/logs/')
        # shell('find {self.output} -not -name "merqury*" -delete')


if __name__ == '__main__':
    Wrapper(snakemake)