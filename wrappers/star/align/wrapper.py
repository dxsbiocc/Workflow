# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/08/20 20:34:52
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import re
import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)

        # same parameters default value
        params = {
            "--outReadsUnmapped": "None",
            "--twopassMode": "Basic",
            "--outSAMstrandField": "intronMotif",       # include for potential use with StringTie for assembly
            "--outSAMunmapped": "Within",
            "--outFilterMultimapScoreRange": 1,
            "--outFilterMultimapNmax": 20,
            "--outFilterMismatchNmax": 10,
            "--sjdbScore": 2,
            "--genomeLoad": "NoSharedMemory",
            "--limitBAMsortRAM": 0,
            "--outFilterMatchNminOverLread": 0.33,
            "--outFilterScoreMinOverLread": 0.33,
            "--sjdbOverhang": 100,
            "--outSAMattributes": "NH HI NM MD AS XS",
            "--outSAMtype": "BAM SortedByCoordinate",
            "--outSAMheaderHD": "@HD VN:1.4",
            # fusion parameters
            "--chimSegmentMin": 12,                     # ** essential to invoke chimeric read detection & reporting **
            "--chimJunctionOverhangMin": 8,
            "--chimOutJunctionFormat": 1,               # **essential** includes required metadata in Chimeric.junction.out file.
            "--alignSJDBoverhangMin": 1,
            "--alignMatesGapMax": 1000000,              # avoid readthru fusions within 100k
            "--alignIntronMax": 500000,
            "--alignSJstitchMismatchNmax": "5 -1 5 5",  # settings improved certain chimera detections
            "--outSAMattrRGline": "ID:GRPundef",
            "--chimMultimapScoreRange": 3,
            "--chimScoreJunctionNonGTAG": -4,
            "--chimMultimapNmax": 20,
            "--chimNonchimScoreDropMin": 10,
            "--peOverlapNbasesMin": 12,
            "--peOverlapMMp": 0.1,
            "--alignInsertionFlush": "Right",
            "--alignSplicedMateMapLminOverLmate": 0,
            "--alignSplicedMateMapLmin": 30
        }
        default_value = ""
        for k, v in params.items():
            if k not in self.extra:
                default_value += f" {k} {v}"
        self.extra += default_value
        # extra patameters
        if "--outSAMtype BAM SortedByCoordinate" in self.extra:
            self.stdout = "BAM_SortedByCoordinate"
        elif "BAM Unsorted" in self.extra:
            self.stdout = "BAM_Unsorted"
        else:
            self.stdout = "SAM"
        # fastq file
        reads = self.snakemake.input.get("reads")
        fq1 = [r for r in reads if '.clean.R1.fq.gz' in r]
        assert fq1, "input-> fastq1 is a required input parameter"
        
        fq2 = [r for r in reads if '.clean.R2.fq.gz' in r]
        if fq2:
            assert len(fq1) == len(fq2), \
                "input-> equal number of files required for fq1 and fastq2"
        # input string
        self.left = ",".join(fq1)
        self.right = ",".join(fq2) if fq2 is not None else ""

        
        # decomplession
        if fq1[0].endswith(".gz"):
            self.readcmd = "--readFilesCommand zcat"
        elif fq1[0].endswith(".bz2"):
            self.readcmd = "--readFilesCommand bzcat"
        else:
            self.readcmd = ""
        
        # genome index
        self.index = self.snakemake.input.get("index")
        if not self.index:
            self.index = self.snakemake.params.get("idx", "")        

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                "STAR "
                " --runThreadN {self.snakemake.threads}"
                " --genomeDir {self.index}"
                " --readFilesIn {self.left} {self.right}"
                " {self.readcmd}"
                " {self.extra}"
                " --outTmpDir {tmpdir}/STARtmp"
                " --outFileNamePrefix {tmpdir}/"
                " --outStd {self.stdout}"
                " > {self.snakemake.output.aln}"
                " {self.log}"
            )

            if self.snakemake.output.get("reads_per_gene"):
                shell("cat {tmpdir}/ReadsPerGene.out.tab > {self.snakemake.output.reads_per_gene:q}")
            if self.snakemake.output.get("chim_junc"):
                shell("cat {tmpdir}/Chimeric.out.junction > {self.snakemake.output.chim_junc:q}")
            if self.snakemake.output.get("chim_junc_sam"):
                shell("cat {tmpdir}/Chimeric.out.sam > {self.snakemake.output.chim_junc_sam:q}")
            if self.snakemake.output.get("sj"):
                shell("cat {tmpdir}/SJ.out.tab > {self.snakemake.output.sj:q}")
            if self.snakemake.output.get("log"):
                shell("cat {tmpdir}/Log.out > {self.snakemake.output.log:q}")
            if self.snakemake.output.get("log_progress"):
                shell("cat {tmpdir}/Log.progress.out > {self.snakemake.output.log_progress:q}")
            if self.snakemake.output.get("log_final"):
                shell("cat {tmpdir}/Log.final.out > {self.snakemake.output.log_final:q}")


if __name__ == '__main__':
    Wrapper(snakemake)