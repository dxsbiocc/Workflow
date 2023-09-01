# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/31 18:20:59
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import tempfile
from os import path
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
        self.range = self.snakemake.params.get("range", " --qrange 0 1")

        self.track = self.snakemake.input.get("track", "")
        track_col_name = self.snakemake.params.get("track_col_name", "")
        if self.track and track_col_name:
            self.track = f" {self.track}::{track_col_name}"

        self.expected = self.snakemake.input.get("expected", "")
        column = self.snakemake.params.get("column", "")
        if column:
            self.expected = f" {self.expected}::{column}"

        self.resolution = self.snakemake.params.get(
            "resolution", self.snakemake.wildcards.get("resolution", 0)
        )
        if not self.resolution:
            raise ValueError("Please specify resolution either as a wildcard or as a parameter")
        

        self.suffix = ""
        self.fig = self.snakemake.output.get("fig", "")
        if self.fig:
            self.suffix = path.splitext(self.fig)[1][1:]
            self.extra += f" --fig {self.suffix}"


    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "(cooltools saddle"
                " {self.snakemake.input.cooler}::resolutions/{self.resolution} "
                " {self.track} "
                " {self.expected} "
                " {self.extra} "
                " {self.range}"
                " -o {tmpdir}/out"
                ") {self.log}"
            )

            shell("mv {tmpdir}/out.saddledump.npz {self.snakemake.output.saddle}")
            shell("mv {tmpdir}/out.digitized.tsv {self.snakemake.output.digitized_track}")
            if self.fig:
                shell("mv {tmpdir}/out.{self.suffix} {self.snakemake.output.fig}")


if __name__ == "__main__":
    Wrapper(snakemake)