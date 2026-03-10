# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2024/01/17 14:35:31
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import subprocess
from snakemake_wrapper_utils.base import WrapperBase, get_logger


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        self.logger = get_logger(snakemake.rule, filename=str(snakemake.log))
        super().__init__(snakemake)

    def parser(self):
        self.bam = self.snakemake.input.get('bam')
        
        self.dnase = self.snakemake.params.get('dnase')
        self.blacklist = self.snakemake.params.get('blacklist')
        self.promoter = self.snakemake.params.get('promoter')
        self.enhancer = self.snakemake.params.get('enhancer')

    def get_output(self, command):
        frac = 0
        try:
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
        except Exception as e:
            self.logger.error(e)
            exit(1)
        else:
            frac = result.stdout.strip()
            self.logger.info(f"get fraction {frac}")
        return frac
    
    
    def run(self):
        cmd = "bedtools sort -i {} | bedtools merge -i stdin | bedtools intersect -u -a {} -b stdin -ubam | samtools view -c"
        total = self.get_output(f"samtools view -c {self.bam}")
        # n_dnase, n_blacklist, n_promoter, n_enhancer
        result = [total, 0, 0, 0, 0]
        if self.dnase:
            result[1] = self.get_output(cmd.format(self.dnase, self.bam))
        if self.blacklist:
            result[2] = self.get_output(cmd.format(self.blacklist, self.bam))
        if self.promoter:
            result[3] = self.get_output(cmd.format(self.promoter, self.bam))
        if self.enhancer:
            result[4] = self.get_output(cmd.format(self.enhancer, self.bam))
        
        with open(str(self.snakemake.output), 'w') as fout:
            fout.write("total\tn_dnase\tn_blacklist\tn_promoter\tn_enhancer\n")
            fout.write('\t'.join(str(i) for i in result))
        self.logger.info(f'result saved in {str(self.snakemake.output)}')

if __name__ == '__main__':
    Wrapper(snakemake)