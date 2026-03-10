# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/09/21 14:50:04
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
        # mapping mode
        self.mode = self.snakemake.params.get("mode", "")
        assert self.mode in ['atac', 'chip', 'hic']
        
        # Single-end read files or paired-end read files
        fq1 = self.snakemake.input.get('fastq1')
        if isinstance(fq1, list):
            fq1 = ','.join(fq1)
        self.extra += f" -1 {fq1}"
        fq2 = self.snakemake.input.get('fastq2', "")
        if fq2:
            if isinstance(fq2, list):
                fq2 = ','.join(fq2)
            self.extra += f" -2 {fq2}"
        # Cell barcode files.
        if self.snakemake.input.get('barcode'):
            self.extra += f" -b {','.join(self.snakemake.input.barcode)}"
        if self.snakemake.input.get('whitelist'):
            self.extra += f" --barcode-whitelist {self.snakemake.input.whitelist}"

        # Remove PCR duplicates
        dedup = self.snakemake.params.get("dedup", "all")
        if dedup == "bulk":
            self.extra += f" --remove-pcr-duplicates-at-bulk-level"
        elif dedup == "cell":
            self.extra += f" --remove-pcr-duplicates-at-cell-level"
        else:
            self.extra += f" --remove-pcr-duplicates"

        # output options
        outfmt = self.snakemake.params.get("outfmt", "")
        if outfmt:
            assert outfmt in ['BED', 'TagAlign', 'pairs', 'SAM']
            self.extra += f" --{outfmt}"

    def run(self):
        shell(
            "chromap"
            " --preset {self.mode}"
            " -r {self.snakemake.input.ref}"
            " -x {self.snakemake.input.index}"
            " -t {self.snakemake.threads}"
            " {self.extra}"
            " -o {self.snakemake.output[0]}"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)