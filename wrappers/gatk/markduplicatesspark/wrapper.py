# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/08 12:26:33
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
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
        self.java_opts = get_java_opts(self.snakemake)
        self.spark_runner = self.snakemake.params.get("spark_runner", "LOCAL")
        self.spark_master = self.snakemake.params.get(
            "spark_master", "local[{}]".format(self.snakemake.threads)
        )
        self.spark_extra = self.snakemake.params.get("spark_extra", "")
        self.java_opts = get_java_opts(self.snakemake)

        metrics = self.snakemake.output.get("metrics", "")
        self.metrics = f"--metrics-file {metrics}" if metrics else ""

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "gatk --java-options '{self.java_opts}' MarkDuplicatesSpark"
                " --input {self.snakemake.input}"
                " {self.extra}"
                " --tmp-dir {tmpdir}"
                " --output {self.snakemake.output.bam}"
                " {self.metrics}"
                " -- --spark-runner {self.spark_runner} --spark-master {self.spark_master} {self.spark_extra}"
                " {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)