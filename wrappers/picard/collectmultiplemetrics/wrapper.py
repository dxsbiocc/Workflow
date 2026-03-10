# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/21 17:14:32
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import os
import pathlib
import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase
from snakemake_wrapper_utils.java import get_java_opts


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.java_opts = get_java_opts(self.snakemake)
        
        self.exts_to_prog = {
            ".alignment_summary_metrics": "CollectAlignmentSummaryMetrics",
            ".insert_size_metrics": "CollectInsertSizeMetrics",
            ".insert_size_histogram.pdf": "CollectInsertSizeMetrics",
            ".quality_distribution_metrics": "QualityScoreDistribution",
            ".quality_distribution.pdf": "QualityScoreDistribution",
            ".quality_by_cycle_metrics": "MeanQualityByCycle",
            ".quality_by_cycle.pdf": "MeanQualityByCycle",
            ".base_distribution_by_cycle_metrics": "CollectBaseDistributionByCycle",
            ".base_distribution_by_cycle.pdf": "CollectBaseDistributionByCycle",
            ".gc_bias.detail_metrics": "CollectGcBiasMetrics",
            ".gc_bias.summary_metrics": "CollectGcBiasMetrics",
            ".gc_bias.pdf": "CollectGcBiasMetrics",
            ".rna_metrics": "RnaSeqMetrics",
            ".bait_bias_detail_metrics": "CollectSequencingArtifactMetrics",
            ".bait_bias_summary_metrics": "CollectSequencingArtifactMetrics",
            ".error_summary_metrics": "CollectSequencingArtifactMetrics",
            ".pre_adapter_detail_metrics": "CollectSequencingArtifactMetrics",
            ".pre_adapter_summary_metrics": "CollectSequencingArtifactMetrics",
            ".quality_yield_metrics": "CollectQualityYieldMetrics",
        }

        # Select programs to run from output files
        progs = set()
        for file in self.snakemake.output:
            matched = False
            for ext in self.exts_to_prog:
                if file.endswith(ext):
                    progs.add(self.exts_to_prog[ext])
                    matched = True
            if not matched:
                raise ValueError(
                    "Unknown type of metrics file requested, for possible metrics files, see \
                        https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/picard/collectmultiplemetrics.html"
                )

        self.programs = "--PROGRAM null --PROGRAM " + " --PROGRAM ".join(progs)


        # Infer common output prefix
        self.output = os.path.commonprefix(self.snakemake.output)

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "picard CollectMultipleMetrics"
                " {self.java_opts} {self.extra}"
                " --INPUT {self.snakemake.input.bam}"
                " --TMP_DIR {tmpdir}"
                " --OUTPUT {self.output}"
                " --REFERENCE_SEQUENCE {self.snakemake.input.ref}"
                " {self.programs}"
                " {self.log}"
            )


        # Under some circumstances, some picard programs might not produce an output 
        # (https://github.com/snakemake/snakemake-wrappers/issues/357). To avoid snakemake 
        # errors, the output files of those programs are created empty (if they do not exist).
        for ext in [
            ext for ext, prog in self.exts_to_prog.items() if prog in ["CollectInsertSizeMetrics"]
        ]:
            for file in self.snakemake.output:
                if file.endswith(ext) and not pathlib.Path(file).is_file():
                    pathlib.Path(file).touch()
            

if __name__ == '__main__':
    Wrapper(snakemake)