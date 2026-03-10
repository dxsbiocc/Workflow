# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2024/04/24 11:01:49
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file
# ============================================================


import pandas as pd
from functools import reduce
from snakemake import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):
    def __init__(self, snakemake):
        super(Wrapper, self).__init__(snakemake)

    def parser(self):
        self.file_list = self.snakemake.input
        if not isinstance(self.file_list, list):
            raise ValueError("The input must be a list of files!")

        self.column_name = self.snakemake.params.get("column_name")
        if not isinstance(self.column_name, list) or len(self.column_name) != len(
            self.file_list
        ):
            raise ValueError("The column name must be a list of strings!")

        self.usecols = self.snakemake.params.get("usecols")
        self.header = self.snakemake.params.get("header")
        # output file
        self.output = str(self.snakemake.output)

    def run(self):
        data = []
        for i, file in enumerate(self.file_list):
            df = pd.read_csv(
                file,
                sep="\t",
                header=self.header,
                usecols=self.usecols,
                index_col=self.usecols[0],
            )
            df.columns = [self.column_name[i]]
            data.append(df)
        data = reduce(
            lambda x, y: pd.merge(x, y, left_index=True, right_index=True), data
        )
        data.to_csv(self.output, sep="\t", index=True)


if __name__ == "__main__":
    Wrapper(snakemake)
