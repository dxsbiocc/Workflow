# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/21 17:09:41
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase
from snakemake_wrapper_utils.java import get_java_opts


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.java_opts = get_java_opts(self.snakemake)
        
        # strand
        self.strand = self.snakemake.params.get("strand", "NONE")
        # reference
        ref = self.snakemake.input.get("ref", "")
        self.ref = f"--REFERENCE_SEQUENCE {ref}" if ref else ""

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "picard CollectRnaSeqMetrics"
                " {self.java_opts} {self.extra}"
                " {self.ref}"
                " --INPUT {self.snakemake.input.bam}"
                " --REF_FLAT {self.snakemake.input.refflat}"
                " --STRAND_SPECIFICITY {self.strand}"
                " --TMP_DIR {tmpdir}"
                " --OUTPUT {self.snakemake.output}"
                " {self.log}"
            )
            

if __name__ == '__main__':
    Wrapper(snakemake)