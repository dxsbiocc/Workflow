# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2024/01/11 20:01:59
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file
# ============================================================


import os
from tempfile import TemporaryDirectory
from snakemake import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.pcl = str(self.snakemake.input.pcl)
        self.top = self.snakemake.params.get("top", 10)

        self.outdir = str(self.snakemake.output)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def run(self):
        # run humann_barplot
        pathway = []
        try:
            with open(self.pcl, "r") as fp:
                text = fp.read().split("\n")
                for line in text[5:]:
                    pathway.append(line.split("\t")[0].split(":")[0])
        except Exception as e:
            print(f"Error reading {self.pcl}: {e}", file=self.log)
            return

        pathway = list(dict.fromkeys(pathway))
        for p in pathway[: self.top]:
            output = os.path.join(self.outdir, f"{p}.pdf")
            shell(
                "humann_barplot"
                " --focal-feature {p}"
                " --input {self.pcl}"
                " --output {output}"
                " {self.extra}"
            )


if __name__ == "__main__":
    Wrapper(snakemake)
