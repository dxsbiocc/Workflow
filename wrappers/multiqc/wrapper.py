# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/28 17:04:31
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
        # Set this to False if multiqc should use the actual input directly
        # instead of parsing the folders where the provided files are located
        use_input_files_only = self.snakemake.params.get("use_input_files_only", False)

        if not use_input_files_only:
            self.input_data = set(os.path.dirname(fp) for fp in self.snakemake.input)
        else:
            self.input_data = set(self.snakemake.input)

        self.output_dir = os.path.dirname(self.snakemake.output[0])
        self.output_name = os.path.basename(self.snakemake.output[0])


    def run(self):
        shell(
            "multiqc"
            " {self.extra}"
            " --force"
            " -o {self.output_dir}"
            " -n {self.output_name}"
            " {self.input_data}"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)