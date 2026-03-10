# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/08 13:47:38
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
        self.annotation_resources = [
            "--resource:{}".format(self.fmt_res(resname, resparams))
            for resname, resparams in self.snakemake.params["resources"].items()
        ]

        self.annotation = list(map("-an {}".format, self.snakemake.params.annotation))
        tranches = self.snakemake.output.get("tranches", "")
        self.tranches = f"--tranches-file {tranches}" if tranches else ""
            
    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "gatk --java-options '{self.java_opts}' VariantRecalibrator"
                " --variant {self.snakemake.input.vcf}"
                " --reference {self.snakemake.input.ref}"
                " --mode {self.snakemake.params.mode}"
                " {self.annotation_resources}"
                " {self.tranches}"
                " {self.annotation}"
                " {self.extra}"
                " --tmp-dir {tmpdir}"
                " --output {self.snakemake.output.vcf}"
                " {self.log}"
            )

    def fmt_res(self, resname, resparams):
        fmt_bool = lambda b: str(b).lower()
        try:
            f = self.snakemake.input.get(resname)
        except KeyError:
            raise RuntimeError(
                f"There must be a named input file for every resource (missing: {resname})"
            )
        return "{},known={},training={},truth={},prior={} {}".format(
            resname,
            fmt_bool(resparams["known"]),
            fmt_bool(resparams["training"]),
            fmt_bool(resparams["truth"]),
            resparams["prior"],
            f,
        )


if __name__ == '__main__':
    Wrapper(snakemake)