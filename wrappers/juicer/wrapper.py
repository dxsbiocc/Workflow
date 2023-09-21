# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/09/06 14:47:36
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
        # output
        self.outdir = os.path.abspath(self.snakemake.output[0])
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        # input
        self.genome = self.snakemake.input.get('genome')
        self.restriction_sites = self.snakemake.input.get('restriction_sites')
        self.chromsizes = self.snakemake.input.get('chromsizes')

        # self.fastq1 = ','.join(self.snakemake.input.get("fastq1"))
        # self.fastq2 = ','.join(self.snakemake.input.get("fastq2"))
        self.fastq = self.snakemake.input.get("fastq")
        if not os.path.exists(os.path.join(self.outdir, 'fastq')):
            os.makedirs(os.path.join(self.outdir, 'fastq'))
        for file in self.fastq:
            basename = os.path.basename(file)
            basename = basename.replace('.R1', '_R1').replace('.R2', '_R2').replace('fq', 'fastq')
            os.symlink(os.path.abspath(file), os.path.join(self.outdir, 'fastq', basename))
        # juicer
        self.juicer = self.snakemake.input.get("juicer", None)
        if not self.juicer:
            self.build_juicer()
            self.juicer = os.path.join(self.outdir, "juicer/")

        self.script = os.path.join(self.juicer, "CPU/juicer.sh")

        # params
        gname = self.snakemake.params.get('gname', 'genome')
        if gname:
            self.extra += f" -g {gname}"
        
        restriction_type = self.snakemake.params.get('restriction_type', '')
        if restriction_type:
            self.extra += f" -s {restriction_type}"

    def build_juicer(self):
        juicer_tools = "https://github.com/aidenlab/Juicebox/releases/download/v2.20.00/juicer_tools.2.20.00.jar"
        juicer = "https://github.com/dxsbiocc/juicer.git"
        command = f"""
        exec > {self.outdir}/build.log
        git clone {juicer}
        mv juicer {self.outdir}/
        wget {juicer_tools} -O juicer_tools.jar
        mv juicer_tools.jar {self.outdir}/juicer/CPU/common
        ln -s {self.outdir}/juicer/CPU {self.outdir}/juicer/scripts
        """
        try:
            shell(command)
        except Exception as e:
            print(e)
            raise Exception("juicer build failed")
    
    def run(self):
        shell(
            "({self.script}"
            " {self.extra}"
            " -t {self.snakemake.threads}"
            " -z {self.genome}"
            " -y {self.restriction_sites}"
            " -p {self.chromsizes}"
            # " -1 {self.fastq1}"
            # " -2 {self.fastq2}"
            " -d {self.outdir}"
            " -D {self.juicer}"
            ") {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)