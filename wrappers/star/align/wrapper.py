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


import os
from tempfile import TemporaryDirectory
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)

        arriba_params = {
            "--outBAMcompression": 0,
            "--peOverlapNbasesMin": 10,
            "--alignSplicedMateMapLminOverLmate": 0.5,
            "--alignSJstitchMismatchNmax": "5 -1 5 5",
            "--chimScoreDropMax": 30,
            "--chimScoreJunctionNonGTAG": 0,
            "--chimScoreSeparation": 1,
            "--chimSegmentReadGapMax": 3,
            "--chimMultimapNmax": 50
        }
        star_fusion_params = {
            "--chimSegmentMin": 12, # ** essential to invoke chimeric read detection & reporting **
            "--chimJunctionOverhangMin": 8,
            "--alignSJstitchMismatchNmax": "5 -1 5 5",
            "--outSAMattrRGline": "ID:GRPundef",
            "--chimMultimapScoreRange": 3,
            "--chimScoreJunctionNonGTAG": -4,
            "--chimMultimapNmax": 20,
            "--chimNonchimScoreDropMin": 10,
            "--peOverlapNbasesMin": 12,
            "--peOverlapMMp": 0.1,
            "--alignInsertionFlush": "Right",
            "--alignSplicedMateMapLminOverLmate": 0,
            "--alignSplicedMateMapLmin": 30,
        }
        tcga_params = {
            "--chimOutJunctionFormat": 1,   # **essential** includes required metadata in Chimeric.junction.out file.
            "--sjdbOverhang": 100,
            "--alignIntronMax": 1000000,
            "--alignIntronMin": 20,
            "--alignMatesGapMax": 1000000,
            "--alignSJDBoverhangMin": 1,
            "--alignSJoverhangMin": 8,
            "--alignSoftClipAtReferenceEnds": "Yes",
            "--chimJunctionOverhangMin": 15,
            "--chimMainSegmentMultNmax": 1,
            "--chimOutType": "Junctions SeparateSAMold WithinBAM SoftClip",
            "--chimSegmentMin": 15,
            "--genomeLoad": "NoSharedMemory",
            "--limitSjdbInsertNsj": 1200000,
            "--outFilterIntronMotifs": "None",
            "--outFilterMatchNminOverLread": 0.33,
            "--outFilterMismatchNmax": 999,
            "--outFilterMismatchNoverLmax": 0.1,
            "--outFilterMultimapNmax": 20,
            "--outFilterScoreMinOverLread": 0.33,
            "--outFilterType": "BySJout",
            "--outSAMattributes": "NH HI NM MD AS XS",
            "--outSAMstrandField": "intronMotif",       # include for potential use with StringTie for assembly
            "--outSAMtype": "BAM SortedByCoordinate",
            "--outSAMunmapped": "Within",
            "--quantMode": "TranscriptomeSAM GeneCounts",
            "--twopassMode": "Basic",
            "--outSAMheaderHD": "@HD VN:1.4",
        }
        
        fusion = self.snakemake.params.get('fusion')
        if fusion and fusion == "arriba":
            tcga_params.update(arriba_params)
        elif fusion and fusion == "star-fusion":
            tcga_params.update(star_fusion_params)
        default_value = ""
        for k, v in tcga_params.items():
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
        fq1 = sorted([r for r in reads if 'R1' in r])
        assert fq1, "input-> fastq1 is a required input parameter"
        
        fq2 = sorted([r for r in reads if 'R2' in r])
        if fq2:
            assert len(fq1) == len(fq2), \
                "input-> equal number of files required for fq1 and fastq2"
        # input string
        self.left = ",".join(fq1)
        self.right = ",".join(fq2) if fq2 is not None else ""

        
        # decomplession
        if fq1[0].endswith(".gz"):
            self.readcmd = "--readFilesCommand gunzip -c"
        elif fq1[0].endswith(".bz2"):
            self.readcmd = "--readFilesCommand bunzip2 -c"
        else:
            self.readcmd = ""

        out_unmapped = self.snakemake.output.get("unmapped", "")
        self.out_unmapped = "--outReadsUnmapped Fastx" if out_unmapped else "--outReadsUnmapped None",
        # genome index
        self.index = self.snakemake.input.get("index")
        if not self.index:
            self.index = self.snakemake.params.get("idx", "")

        self.output = os.path.dirname(self.snakemake.output.aln)

    def run(self):
        with TemporaryDirectory() as tmpdir:
            shell(
                "STAR "
                " --runThreadN {self.snakemake.threads}"
                " --genomeDir {self.index}"
                " --readFilesIn {self.left} {self.right}"
                " {self.readcmd}"
                " {self.out_unmapped}"
                " {self.extra}"
                " --outTmpDir {tmpdir}/STARtmp"
                " --outFileNamePrefix {tmpdir}/"
                " --outStd {self.stdout}"
                " > {self.snakemake.output.aln}"
                " {self.log}"
            )
            # unmapped reads
            unmapped = list(self.snakemake.output.get("unmapped"))
            if unmapped:
                for i, out_unmapped in enumerate(unmapped, 1):
                    if os.path.exists(f"{tmpdir}/Unmapped.out.mate{i}"):
                        cmd = "gzip -c" if out_unmapped.endswith("gz") else "cat"
                        shell("{cmd} {tmpdir}/Unmapped.out.mate{i} > {out_unmapped}")
                    else:
                        shell("touch {out_unmapped}")
            if not os.path.exists(self.output):
                os.makedirs(self.output)
            # There have format string in shell function, so it can't be used
            # shell(f"find {tmpdir} -maxdepth 1 -type f -exec mv {{}} {self.output} \\;")
            cmd = f"find {tmpdir} -maxdepth 1 -type f -exec mv {{}} {self.output} \\;"
            os.system(cmd)


if __name__ == '__main__':
    Wrapper(snakemake)