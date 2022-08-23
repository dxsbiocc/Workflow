# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 10:27:29
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        # java options
        self.java_opts = get_java_opts(self.snakemake)
        self.log = self.snakemake.log_fmt_shell()

        bams = self.snakemake.input.bams
        if isinstance(bams, str):
            bams = [bams]
        self.bams = list(map("--INPUT {}".format, bams))

        if self.snakemake.output.bam.endswith(".cram"):
            self.output = "/dev/stdout"
            if self.snakemake.params.embed_ref:
                view_options = "-O cram,embed_ref"
            else:
                view_options = "-O cram"
            self.convert = f" | samtools view -@ {self.snakemake.threads} {view_options} \
                            --reference {self.snakemake.input.ref} -o {self.snakemake.output.bam}"
        else:
            self.output = self.snakemake.output.bam
            self.convert = ""

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "(picard MarkDuplicates"  # Tool and its subcommand
                " {self.java_opts}"  # Automatic java option
                " {self.extra}"  # User defined parmeters
                " {self.bams}"  # Input bam(s)
                " --TMP_DIR {tmpdir}"
                " --OUTPUT {self.output}"  # Output bam
                " --METRICS_FILE {self.snakemake.output.metrics}"  # Output metrics
                " {self.convert} ) {self.log}"  # Logging
            )


if __name__ == '__main__':
    Wrapper(snakemake)