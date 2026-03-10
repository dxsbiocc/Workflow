# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 09:29:47
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file
# ============================================================


import tempfile
from pathlib import Path
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase
from snakemake_wrapper_utils.snakemake import get_mem


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.memory_requirements = get_mem(self.snakemake, out_unit="KiB") * 1024
        self.log = self.snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

        # parse short reads
        reads = (
            self.snakemake.input.reads
            if hasattr(self.snakemake.input, "reads")
            else self.snakemake.input
        )

        input_arg = ""

        # handle named inputs if available
        if hasattr(self.snakemake.input, "r1") and hasattr(self.snakemake.input, "r2"):
            if isinstance(self.snakemake.input.r1, str):
                self.snakemake.input.r1 = [self.snakemake.input.r1]
            if isinstance(self.snakemake.input.r2, str):
                self.snakemake.input.r2 = [self.snakemake.input.r2]
            input_arg += f" -1 {','.join(self.snakemake.input.r1)} -2 {','.join(self.snakemake.input.r2)}"
        elif len(reads) >= 2:
            input_arg += f" -1 {reads[0]} -2 {reads[1]}"

        # handle interleaved reads if specified
        if hasattr(self.snakemake.input, "interleaved"):
            input_arg += f" --12 {self.snakemake.input.interleaved}"
        elif len(reads) >= 3 and not hasattr(self.snakemake.input, "r1"):
            input_arg += f" --12 {reads[2]}"

        # handle additional reads if specified
        if hasattr(self.snakemake.input, "unpaired"):
            input_arg += f" --read {self.snakemake.input.unpaired}"
        elif len(reads) >= 4 and not hasattr(self.snakemake.input, "r1"):
            input_arg += f" --read {reads[3]}"

        self.input_arg = input_arg

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_tmpdir = Path(tmpdir) / "temp"

            shell(
                "megahit"
                " -t {self.snakemake.threads}"
                " -m {self.memory_requirements}"
                " -o {output_tmpdir}"
                " {self.input_arg}"
                " {self.extra}"
            )

            # Ensure user can name each file according to their need
            output_mapping = {
                "contigs": f"{output_tmpdir}/final.contigs.fa",
                "json": f"{output_tmpdir}/options.json",
                "log": f"{output_tmpdir}/log",
            }
            for output_key, temp_file in output_mapping.items():
                output_path = self.snakemake.output.get(output_key)
                if output_path:
                    shell("cp --verbose {temp_file:q} {output_path:q} {self.log}")


if __name__ == "__main__":
    Wrapper(snakemake)
