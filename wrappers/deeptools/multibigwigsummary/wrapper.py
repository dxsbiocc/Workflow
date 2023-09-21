# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/10/20 17:24:26
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
        self.subcommand = self.snakemake.params.get('subcommand')
        assert self.subcommand in ["bins", "BED-file"], \
            f"subcommand '{self.subcommand}' must be 'bins' or 'BED-file'"
        self.optional = ''
        if self.subcommand == 'BED-file':
            self.optional += f'--BED {self.snakemake.input.bed} '
        
        raw_count = self.snakemake.output.get('raw_count')
        self.optional += f'--outRawCounts {raw_count} ' if raw_count else ''

        labels = self.snakemake.params.get("labels", "")
        if labels:
            labels = ' '.join([os.path.basename(label).replace(' ', '-') for label in labels])
            self.optional += f" --labels {labels} "

    def run(self):
        shell(
            "( multiBigwigSummary "
            "{self.subcommand} "
            "--bwfiles {self.snakemake.input.bigwig} "
            "--outFileName {self.snakemake.output.out} "
            "--numberOfProcessors {self.snakemake.threads} "
            "{self.optional} "
            "{self.extra} "
            ") "
            "{self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)