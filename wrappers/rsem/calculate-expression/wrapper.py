# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/04/21 14:12:19
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import os
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        # input files
        bam = self.snakemake.input.get("bam", "")
        fq_one = self.snakemake.input.get("fq_one", "")
        fq_two = self.snakemake.input.get("fq_two", "")
        # mode
        mode = self.snakemake.params.get("mode", "fastq")
        assert mode in ["fastq", "bam"], "You must choose one of 'bam' or 'fastq'"
        if mode == "bam":
            self.input_bam = "--alignments"
            self.input_string = bam
            paired_end = self.snakemake.params.get("paired_end", False)
        else:
            self.input_bam = ""
            if fq_one:
                if isinstance(fq_one, list):
                    num_fq_one = len(fq_one)
                    self.input_string = ",".join(fq_one)
                    if all([1 if fq.endswith('gz') else 0 for fq in fq_one]):
                        self.extra += " --star-gzipped-read-file"
                    elif all([1 if fq.endswith('bz2') else 0 for fq in fq_one]):
                        self.extra += " --star-bzipped-read-file"
                else:
                    num_fq_one = 1
                    self.input_string = fq_one
                    if fq_one.endswith('gz'):
                        self.extra += " --star-gzipped-read-file"
                    elif fq_one.endswith('bz2'):
                        self.extra += " --star-bzipped-read-file"
                
                if fq_two:
                    paired_end = True
                    if isinstance(fq_two, list):
                        num_fq_two = len(fq_two)
                        if num_fq_one != num_fq_two:
                            raise Exception(
                                "Got {} R1 FASTQs, {} R2 FASTQs.".format(num_fq_one, num_fq_two)
                            )
                    else:
                        fq_two = [fq_two]
                    self.input_string += " " + ",".join(fq_two)
                else:
                    paired_end = False
        if paired_end:
            self.paired_end_string = "--paired-end"
        else:
            self.paired_end_string = ""

        self.output = self.snakemake.output.get('quant')
        outdir = os.path.dirname(self.output)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        self.outprefix = os.path.join(outdir, self.snakemake.wildcards.sample)
        
        mapping = self.snakemake.params.get("mapping")
        if bam and mapping:
            mapping = False
        if mapping == "bowtie2":
            self.extra += " --bowtie2"
        elif mapping == "star":
            self.extra += " --star"
        elif mapping == "hisat2":
            self.extra += " --hisat2-hca"

        self.reference_prefix = os.path.splitext(self.snakemake.input.reference[0])[0]

    def run(self):
        shell(
            "rsem-calculate-expression"
            " --num-threads {self.snakemake.threads}"
            " {self.extra}"
            " {self.paired_end_string}"
            " {self.input_bam}"
            " {self.input_string}"
            " {self.reference_prefix}"
            " {self.outprefix}"
            " {self.log}"
        )
        shell("mv {self.outprefix}.genes.results {self.output}")


if __name__ == '__main__':
    Wrapper(snakemake)