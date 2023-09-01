# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/08/31 20:41:02
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):
    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        binsize = self.snakemake.params.get("binsize", None)
        if binsize:
            self.extra += f" --binSize {' '.join(map(str, binsize))}"
        region = self.snakemake.params.get("region", None)
        if region:
            self.extra += f" --region {region}"
        chromosomeSizes = self.snakemake.params.get("chromosomeSizes", None)
        if chromosomeSizes:
            self.extra += f" --chromosomeSizes {chromosomeSizes}"

        genomeAssembly = self.snakemake.params.get("genomeAssembly", None)
        if genomeAssembly:
            self.extra += f" --genomeAssembly {genomeAssembly}"
    
    def run(self):
        shell(
            "(hicBuildMatrix"
            " -s {self.snakemake.input.bam}"
            " -o {self.snakemake.output.matrix}"
            " --QCfolder {self.snakemake.params.qc}"
            " --restrictionCutFile {self.snakemake.input.rs}"
            " --restrictionSequence {self.snakemake.params.seq}"
            " --danglingSequence {self.snakemake.params.dseq}"
            " --threads {self.snakemake.threads}"
            " {self.extra}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)