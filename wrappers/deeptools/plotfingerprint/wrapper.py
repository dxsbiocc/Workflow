# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/19 15:33:51
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import re
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):

        jsd_sample = self.snakemake.input.get("jsd_sample", None)
        self.jsd = " --JSDsample {jsd} ".format(jsd=jsd_sample) if jsd_sample else ""
        
        out_counts = self.snakemake.output.get("counts")
        self.out_metrics = self.snakemake.output.get("qc_metrics")
        
        self.optional_output = ""
        if out_counts:
            self.optional_output += " --outRawCounts {out_counts} ".format(out_counts=out_counts)

        if self.out_metrics:
            self.optional_output += " --outQualityMetrics {metrics} ".format(metrics=self.out_metrics)

    def run(self):
        shell(
            "(plotFingerprint "
            "-b {self.snakemake.input.bams} "
            "-o {self.snakemake.output.fingerprint} "
            "{self.optional_output} "
            "--numberOfProcessors {self.snakemake.threads} "
            "{self.jsd} "
            "{self.extra}) {self.log}"
        )
        # ToDo: remove the 'NA' string replacement when fixed in deepTools, see:
        # https://github.com/deeptools/deepTools/pull/999
        regex_passes = 2

        with open(self.out_metrics, "rt") as f:
            metrics = f.read()
            for i in range(regex_passes):
                metrics = re.sub("\tNA(\t|\n)", "\tnan\\1", metrics)

        with open(self.out_metrics, "wt") as f:
            f.write(metrics)


if __name__ == "__main__":
    Wrapper(snakemake)