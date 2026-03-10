# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2024/01/18 17:24:52
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from  snakemake import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)

        self.bam = str(self.snakemake.input)
        self.outfile = str(self.snakemake.output)
        self.paired = self.snakemake.params.get('paired')
        self.mito_name = self.snakemake.params.get('mito_name', 'chrM')
    
    def run(self):
        if self.paired:
            shell(
                'echo -e "Total Fragments\tDistinct Fragments\tPositions with One Read\tPositions with Two Read\tNRF = Distinct/Total\tPBC1 = OneRead/Distinct\tPBC2 = OneRead/TwoRead" > {self.outfile};'
                'samtools sort -n -@ 4 {self.bam} -o - |'
                'bedtools bamtobed -bedpe -i - | '
                'awk \'BEGIN{{OFS="\\t"}}{{print $1,$2,$4,$6,$9,$10}}\' | '
                'grep -v "^{self.mito_name}\\s" | sort | uniq -c | '
                'awk \'BEGIN{{mt=0;m0=0;m1=0;m2=0}} ($1==1){{m1=m1+1}} '
                '($1==2){{m2=m2+1}} {{m0=m0+1}} '
                '{{mt=mt+$1}} END{{m1_m2=-1.0; '
                'if(m2>0) m1_m2=m1/m2; m0_mt=0; '
                'if (mt>0) m0_mt=m0/mt; m1_m0=0; if (m0>0) m1_m0=m1/m0; '
                'printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n"'
                ',mt,m0,m1,m2,m0_mt,m1_m0,m1_m2}}\' >> {self.outfile} '
                '{self.log}'
            )
        else:
            shell(
                'echo -e "Total Fragments\tDistinct Fragments\tPositions with One Read\tPositions with Two Read\tNRF = Distinct/Total\tPBC1 = OneRead/Distinct\tPBC2 = OneRead/TwoRead" > {self.outfile};'
                'samtools sort -n -@ 4 {self.bam} -o - |'
                'bedtools bamtobed -bedpe -i - | '
                'awk \'BEGIN{{OFS="\\t"}}{{print $1,$2,$3,$6}}\' | '
                'grep -v "^{self.mito_name}\\s" | sort | uniq -c | '
                'awk \'BEGIN{{mt=0;m0=0;m1=0;m2=0}} ($1==1){{m1=m1+1}} '
                '($1==2){{m2=m2+1}} {{m0=m0+1}} '
                '{{mt=mt+$1}} END{{m1_m2=-1.0; '
                'if(m2>0) m1_m2=m1/m2; m0_mt=0; '
                'if (mt>0) m0_mt=m0/mt; m1_m0=0; if (m0>0) m1_m0=m1/m0; '
                'printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n"'
                ',mt,m0,m1,m2,m0_mt,m1_m0,m1_m2}}\' >> {self.outfile} '
                '{self.log}'
            )
    

if __name__ == '__main__':
    Wrapper(snakemake)