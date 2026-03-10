# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2024/04/24 11:01:49
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
    def __init__(self, snakemake):
        super(Wrapper, self).__init__(snakemake)

    def parser(self):
        self.bam = self.snakemake.input.get('bam')
        self.bed = self.snakemake.input.get('bed')
        self.outfile = str(self.snakemake.output)
    
    def run(self):
        cmd = f"""
        set -euo pipefail;
        echo -e "Total reads\tPeak reads\tRatio" > {self.outfile};
        # calculate FRiP
        total=$(samtools view -c {self.bam})
        if [ "$total" -eq 0 ]; then
            echo -e "No reads found in the BAM file." > "{self.snakemake.log}"
            exit 1
        fi
        pc=$(bedtools sort -i {self.bed} | bedtools merge -i stdin | bedtools intersect -u -a {self.bam} -b stdin -ubam | samtools view -c)
        ratio=$(echo "scale=4; $pc / $total" | bc)
        # write output
        echo -e "$total\t$pc\t$ratio" >> {self.outfile}
        """
        shell(cmd)


if __name__ == '__main__':
    Wrapper(snakemake)