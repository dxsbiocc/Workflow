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


from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.optional = ''
        corr = self.snakemake.params.get('cor', 'spearman')
        assert corr in ['pearson', 'spearman'], f"Correlation method '{corr}' not support, possible choices: spearman, pearson"
        self.optional += f"--corMethod {corr} "
        
        title = self.snakemake.params.get('title', '')
        if 'spearman' in title and corr != 'title':
            title = title.replace('Spearman', corr)
        self.optional += f'--plotTitle "{title}" '
        
        labels = self.snakemake.params.get('labels')
        self.optional += f"--labels {' '.join(labels)} " if labels else ""

        mat = self.snakemake.output.get('matrix')
        self.optional += f"--outFileCorMatrix {mat} " if mat else ""


    def run(self):
        shell(
            "( plotCorrelation "
            "--corData {self.snakemake.input[0]} "
            "--plotFile {self.snakemake.output.img} "
            "{self.optional} "
            "{self.extra} "
            ") "
            "{self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)