# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/08 15:28:40
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import os
import tempfile
from pathlib import Path
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.java_opts = get_java_opts(self.snakemake)

        self.n_out_files = len(self.snakemake.output)
        assert self.n_out_files > 1, "you need to specify more than 2 output files!"

        self.prefix = Path(os.path.commonprefix(self.snakemake.output))
        self.suffix = os.path.commonprefix([file[::-1] for file in self.snakemake.output])[::-1]
        chunk_labels = [
            out.removeprefix(str(self.prefix)).removesuffix(self.suffix) for out in self.snakemake.output
        ]
        assert all(
            [chunk_label.isnumeric() for chunk_label in chunk_labels]
        ), "all chunk labels have to be numeric!"
        
        self.len_chunk_labels = set([len(chunk_label) for chunk_label in chunk_labels])
        assert len(self.len_chunk_labels) == 1, "all chunk labels must have the same length!"

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "gatk --java-options '{self.java_opts}' SplitIntervals"
                " --intervals {self.snakemake.input.intervals}"
                " --reference {self.snakemake.input.ref}"
                " --scatter-count {self.n_out_files}"
                " {self.extra}"
                " --tmp-dir {tmpdir}"
                " --output {self.prefix.parent}"
                " --interval-file-prefix {self.prefix.name:q}"
                " --interval-file-num-digits {self.len_chunk_labels}"
                " --extension {self.suffix:q}"
                " {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)