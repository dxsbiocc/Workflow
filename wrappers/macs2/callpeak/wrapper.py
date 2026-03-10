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


import os
import sys
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        # control input
        in_contr = self.snakemake.input.get("control", None)
        self.opt_input = "-c {contr}".format(contr=in_contr) if in_contr else ""
        # output
        # ext = "_peaks.xls"
        out_file = os.path.commonprefix(self.snakemake.output)[:-1]
        self.out_name = f'-n {os.path.basename(out_file)}' if out_file else ""
        self.out_dir = f"--outdir {os.path.dirname(out_file)}" if out_file else ""
            
        # verify parameters
        if any(out.endswith(("_peaks.narrowPeak", "_summits.bed")) for out in self.snakemake.output):
            if any(
                out.endswith(("_peaks.broadPeak", "_peaks.gappedPeak"))
                for out in self.snakemake.output
            ):
                sys.exit(
                    "Output files with _peaks.narrowPeak and/or _summits.bed extensions cannot be created together with _peaks.broadPeak and/or _peaks.gappedPeak extended output files.\n"
                    "For usable extensions please see https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.\n"
                )
            else:
                if " --broad" in self.extra:
                    sys.exit(
                        "If --broad option in params is given, the _peaks.narrowPeak and _summits.bed files will not be created. \n"
                        "Remove --broad option from params if these files are needed.\n"
                    )
        if any(
            out.endswith(("_peaks.broadPeak", "_peaks.gappedPeak")) for out in self.snakemake.output
        ):
            # exclude --broad-cutoff parameter
            if "--broad " not in self.extra and not self.extra.endswith("--broad"):
                self.extra += " --broad "
        if any(
            out.endswith(("_treat_pileup.bdg", "_control_lambda.bdg"))
            for out in self.snakemake.output
        ):
            if all(p not in self.extra for p in ["--bdg", "-B"]):
                self.extra += " --bdg "
        else:
            if any(p in self.extra for p in ["--bdg", "-B"]):
                sys.exit(
                    "If --bdg or -B option in params is given, the _control_lambda.bdg and _treat_pileup.bdg extended files must be specified in output. \n"
                )

    def run(self):
        shell(
            "(macs2 callpeak "
            "-t {self.snakemake.input.treatment} "
            "{self.opt_input} "
            "{self.out_dir} "
            "{self.out_name} "
            "{self.extra}) {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)