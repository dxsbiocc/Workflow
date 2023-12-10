# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/31 16:54:52
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)
        ## Extract arguments
        view = self.snakemake.input.get("view", "")
        if view:
            self.extra += f" --view {view}"
        track = self.snakemake.input.get("track", "")
        track_col_name = self.snakemake.params.get("track_col_name", "")
        if track and track_col_name:
            self.extra += f" --phasing-track {track}::{track_col_name}"
        elif track:
            self.extra += f" --phasing-track {track}"


        self.bigwig = "--bigwig" if self.snakemake.output.get("bigwig", "") else ""

        self.resolution = self.snakemake.params.get("resolution", self.snakemake.wildcards.get("resolution"))
        assert (
            self.resolution
        ), "Please specify resolution either as a `wildcard` or as a `parameter`"

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "cooltools eigs-trans"
                " {self.extra} "
                " {self.bigwig}"
                " -o {tmpdir}/out"
                " {self.snakemake.input.cooler}::resolutions/{self.resolution} "
                " {self.log}"
            )

            shell("mv {tmpdir}/out.trans.vecs.tsv {self.snakemake.output.vecs}")
            shell("mv {tmpdir}/out.trans.lam.txt {self.snakemake.output.lam}")
            if self.bigwig:
                shell("mv {tmpdir}/out.trans.bw {self.snakemake.output.bigwig}")



if __name__ == "__main__":
    Wrapper(snakemake)
