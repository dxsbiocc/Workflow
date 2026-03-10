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
import sys
import subprocess
from tempfile import TemporaryDirectory
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True)

        # unmapped reads
        out_unmapped = self.snakemake.output.get("unmapped", "")
        self.out_unmapped = "--outReadsUnmapped Fastx" if out_unmapped else "--outReadsUnmapped None",
        # genome index
        self.index = self.snakemake.input.get("index")
        if not self.index:
            self.index = self.snakemake.params.get("idx", "")

        self.output = os.path.abspath(str(self.snakemake.output))
        # fastq file
        self.fq1 = self.snakemake.input.get("fq1", "")
        self.fq2 = self.snakemake.input.get("fq2", "")
        file = self.fq1 if isinstance(self.fq1, str) else self.fq1[0]
        if self.fq1:
            if isinstance(self.fq1, list):
                self.fq1 = ','.join(self.fq1)
                if self.fq2:
                    self.fq2 = ','.join(self.fq2)

        self.input = f"--readFilesIn {self.fq2} {self.fq1}"
        # decomplession
        if file.endswith(".gz"):
            self.readcmd = "--readFilesCommand gunzip -c"
            self.zcmd = "zcat"
        elif file.endswith(".bz2"):
            self.readcmd = "--readFilesCommand bunzip2 -c"
            self.zcmd = "bzcat"
        else:
            self.readcmd = ""
            self.zcmd = "cat"
        # params
        self.barcode = self.snakemake.params.get("barcode")
        cblen = self.snakemake.params.get("cblen", 8)
        umilen = self.snakemake.params.get("umilen", 8)

        params = {
            "runDirPerm": "All_RWX",
            "outSAMtype": "BAM Unsorted",
            'soloFeatures': 'Gene GeneFull',
            # 'outSAMattributes': 'NH HI AS nM CB UB CR CY UR UY GX GN',  # only be output in the sorted BAM file
            # 'outBAMsortingBinsN': 500,
            # 'limitBAMsortRAM': 60000000000,
            # 'outMultimapperOrder': 'Random',
            # 'runRNGseed': 1,
        }

        self.scrna = self.snakemake.params.get("scrna", False)
        if self.scrna == '10x':
            # auto infer barcode and umi length
            self.auto_10x()
            # update params
            params.update({
                'soloType': 'CB_UMI_Simple',
                'soloCBwhitelist': os.path.join(self.barcode, self.scrna, self.BC),
                'soloCBlen': self.cblen,
                'soloUMIstart': self.cblen + 1,
                'soloUMIlen': self.umilen,
                'soloUMIdedup': '1MM_CR',
                'soloCBmatchWLtype': '1MM_multi_Nbase_pseudocounts',
                'soloUMIfiltering': 'MultiGeneUMI_CR',
                'soloCellFilter': 'EmptyDrops_CR',
                'outFilterScoreMin': 30,
                'soloFeatures': 'Gene GeneFull Velocyto',
                'soloMultiMappers': 'EM',
            })
            if self.paired:
                params.update({
                    'soloBarcodeMate': 1,
                    'clip5pNbases': '39 0',
                    'soloCBstart': 1,
                    'soloStrand': 'Forward',
                    'outFilterScoreMin': 30,
                })
                self.input = f"--readFilesIn {self.fq1} {self.fq2}"
            else:
                params.update({
                    'soloBarcodeReadLength': 0,
                    'soloStrand': self.strand,
                    'clipAdapterType': 'CellRanger4',
                })
        elif self.scrna == 'drop':
            params.update({
                'soloType': 'CB_UMI_Simple',
                'soloCBwhitelist': 'None',
                'soloCBstart': 1,
                'soloCBlen': cblen,         # 12
                'soloUMIstart': cblen + 1,
                'soloUMIlen': umilen,       # 8
                'soloBarcodeReadLength': 0,
            })
        elif self.scrna == 'indrop':
            bclist = [os.path.join(self.barcode, self.scrna, f'gel_barcode{i}_list.txt') for i in range(1, 3)]
            adapter = self.snakemake.params.get("adapter", "GAGTGATTGCTTGTGACGCCTT")
            params.update({
                'soloCBwhitelist': ' '.join(bclist),
                'soloAdapterSequence': adapter,
                'soloType': 'CB_UMI_Complex',
                'soloAdapterMismatchesNmax': 3,
                'soloCBmatchWLtype': '1MM',
                'soloCBposition': '0_0_2_-1 3_1_3_8',
                'soloUMIposition': '3_9_3_14',
            })
        elif self.scrna == 'smart2':
            params.update({
                'soloType': 'SmartSeq',
                'soloUMIdedup': 'Exact',
                'soloStrand': 'Unstranded',
                'limitOutSJcollapsed': 10000000,
                'soloCellFilter': 'None',
            })
        elif self.scrna == 'bd_rhapsody':
            bclist = [os.path.join(self.barcode, self.scrna, f'bd_rhapsody_barcode{i}.txt') for i in range(1, 4)]
            params.update({
                'soloCBwhitelist': ' '.join(bclist),
                'soloType': 'CB_UMI_Complex',
                'soloUMIlen': umilen,  # 8
                'soloCBmatchWLtype': '1MM',
                'soloCBposition': '0_0_0_8 0_21_0_29 0_43_0_51',
                'soloUMIposition': '0_52_0_59',
            })
        elif self.scrna == 'cel':
            params.update({
                'soloCBwhitelist': 'None',
                'soloType': 'CB_UMI_Simple',
                'clip5pNbases': '12 0',
                'soloUMIstart': 1,
                'soloUMIlen': umilen,  # 6
                'soloCBstart': umilen + 1,
                'soloCBlen': cblen,    # 6
                'soloUMIdedup': 'Exact',
            })
        else:
            raise ValueError(f"`scrna` must be one of 10x, drop, indrop, smart2, strt, bd_rhapsody, cel")

        default_value = ""
        for k, v in params.items():
            if k not in self.extra:
                default_value += f" --{k} {v}"
        self.extra += default_value

    def run_command(self, cmd):
        """Run a shell command and return the output."""
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            sys.exit(result.returncode)
        return result.stdout

    def auto_10x(self):

        with TemporaryDirectory() as tmpdir:
            for i, elem in enumerate(self.fq1.split(',')):
                shell.check_output(f"{self.zcmd} {elem} | head -4000000 > {tmpdir}/{i}.R1_head")
            for i, elem in enumerate(self.fq2.split(',')):
                shell.check_output(f"{self.zcmd} {elem} | head -4000000 > {tmpdir}/{i}.R2_head")

            shell.check_output(f"cat {tmpdir}/*.R1_head | seqtk sample -s100 - 200000 > {tmpdir}/test.R1.fastq")
            shell.check_output(f"cat {tmpdir}/*.R2_head | seqtk sample -s100 - 200000 > {tmpdir}/test.R2.fastq")
            # rm temp files
            shell(f"rm {tmpdir}/*.R1_head {tmpdir}/*.R2_head")
            # Evaluate the reads for barcode whitelist
            nbc1 = int(self.run_command(f"cat {tmpdir}/test.R1.fastq | awk 'NR%4==2' | cut -c-14 | grep -F -f {self.barcode}/{self.scrna}/737K-april-2014_rc.txt | wc -l").strip())
            nbc2 = int(self.run_command(f"cat {tmpdir}/test.R1.fastq | awk 'NR%4==2' | cut -c-16 | grep -F -f {self.barcode}/{self.scrna}/737K-august-2016.txt | wc -l").strip())
            nbc3 = int(self.run_command(f"cat {tmpdir}/test.R1.fastq | awk 'NR%4==2' | cut -c-16 | grep -F -f {self.barcode}/{self.scrna}/3M-february-2018.txt | wc -l").strip())
            nbcA = int(self.run_command(f"cat {tmpdir}/test.R1.fastq | awk 'NR%4==2' | cut -c-16 | grep -F -f {self.barcode}/{self.scrna}/737K-arc-v1.txt | wc -l").strip())
            R1LEN = int(self.run_command(f"cat {tmpdir}/test.R1.fastq | awk 'NR%4==2' | awk '{{sum+=length($0)}} END {{printf \"%d\\n\",sum/NR+0.5}}'").strip())
            R2LEN = int(self.run_command(f"cat {tmpdir}/test.R2.fastq | awk 'NR%4==2' | awk '{{sum+=length($0)}} END {{printf \"%d\\n\",sum/NR+0.5}}'").strip())
            R1DIS = int(self.run_command(f"cat {tmpdir}/test.R1.fastq | awk 'NR%4==2' | awk '{{print length($0)}}' | sort | uniq -c | wc -l").strip())
            # Choose the right barcode whitelist
            self.BC = ""
            if nbc3 > 50000:
                self.BC = f"3M-february-2018.txt"
            elif nbc2 > 50000:
                self.BC = f"737K-august-2016.txt"
            elif nbcA > 50000:
                self.BC = f"737K-arc-v1.txt"
            elif nbc1 > 50000:
                self.BC = f"737K-april-2014_rc.txt"
            else:
                print("ERROR: No whitelist has matched a random selection of 200,000 barcodes! Match counts: {} (v1), {} (v2), {} (v3), {} (multiome).".format(nbc1, nbc2, nbc3, nbcA), file=sys.stderr)
                sys.exit(1)

            # Check read lengths and set variables
            self.paired = False
            self.umilen, self.cblen = 0, 0
            if R1DIS > 1 and R1LEN <= 30:
                print("ERROR: Read 1 (barcode) has varying length; possibly someone thought it's a good idea to quality-trim it. Please check the fastq files.", file=sys.stderr)
                sys.exit(1)
            elif R1LEN < 24:
                print("ERROR: Read 1 (barcode) is less than 24 bp in length. Please check the fastq files.", file=sys.stderr)
                sys.exit(1)
            elif R2LEN < 40:
                print("ERROR: Read 2 (biological read) is less than 40 bp in length. Please check the fastq files.", file=sys.stderr)
                sys.exit(1)

            if R1LEN > 50:
                self.paired = True

            if self.BC == f"3M-february-2018.txt" or self.BC == f"737K-arc-v1.txt":
                self.cblen, self.umilen = 16, 12
            elif self.BC == f"737K-august-2016.txt":
                self.cblen, self.umilen = 16, 10
            elif self.BC == f"737K-april-2014_rc.txt":
                self.cblen, self.umilen = 14, 10

            if self.cblen + self.umilen > R1LEN:
                NEWUMI = R1LEN - self.cblen
                BCUMI = self.umilen + self.cblen
                print(f"WARNING: Read 1 length ({R1LEN}) is less than the sum of appropriate barcode and UMI ({BCUMI}). Changing UMI setting from {self.umilen} to {NEWUMI}!", file=sys.stderr)
                self.umilen = NEWUMI
            elif self.cblen + self.umilen < R1LEN:
                BCUMI = self.umilen + self.cblen
                print(f"WARNING: Read 1 length ({R1LEN}) is more than the sum of appropriate barcode and UMI ({BCUMI}).", file=sys.stderr)

            # Determine the strand
            self.strand = "Forward"
            shell(f"STAR --runThreadN {self.snakemake.threads} --genomeDir {self.index} --readFilesIn {tmpdir}/test.R2.fastq {tmpdir}/test.R1.fastq --runDirPerm All_RWX --outSAMtype None "
                  f"--soloType CB_UMI_Simple --soloCBwhitelist {self.barcode}/{self.scrna}/{self.BC} --soloBarcodeReadLength 0 --soloCBlen {self.cblen} --soloUMIstart {self.cblen+1} "
                  f"--soloUMIlen {self.umilen} --soloStrand Forward "
                  f"--soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR "
                  f"--soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloFeatures Gene GeneFull "
                  f"--outTmpDir {tmpdir}/STARtmp --outFileNamePrefix {tmpdir} --soloOutFileNames test_forward/ features.tsv barcodes.tsv matrix.mtx {self.log}")

            shell(f"STAR --runThreadN {self.snakemake.threads} --genomeDir {self.index} --readFilesIn {tmpdir}/test.R2.fastq {tmpdir}/test.R1.fastq --runDirPerm All_RWX --outSAMtype None "
                  f"--soloType CB_UMI_Simple --soloCBwhitelist {self.barcode}/{self.scrna}/{self.BC} --soloBarcodeReadLength 0 --soloCBlen {self.cblen} --soloUMIstart {self.cblen+1} "
                  f"--soloUMIlen {self.umilen} --soloStrand Reverse "
                  f"--soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR "
                  f"--soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloFeatures Gene GeneFull "
                  f"--outTmpDir {tmpdir}/STARtmp --outFileNamePrefix {tmpdir} --soloOutFileNames test_reverse/ features.tsv barcodes.tsv matrix.mtx {self.log}")

            fwd = self.run_command(f"grep 'Reads Mapped to GeneFull: Unique GeneFull' {tmpdir}/test_forward/GeneFull/Summary.csv | awk -F ',' '{{printf \"%d\\n\",$2*100+0.5}}'").strip()
            PCTFWD = int(fwd) if fwd.isnumeric() else 0
            rev = self.run_command(f"grep 'Reads Mapped to GeneFull: Unique GeneFull' {tmpdir}/test_reverse/GeneFull/Summary.csv | awk -F ',' '{{printf \"%d\\n\",$2*100+0.5}}'").strip()
            PCTREV = int(rev) if rev.isnumeric() else 0

            if PCTREV >= PCTFWD:
                self.strand = "Reverse"

            if PCTREV < 50 and PCTFWD < 50:
                print(f"WARNING: Low percentage of reads mapping to GeneFull: forward = {PCTFWD} , reverse = {PCTREV}", file=sys.stderr)

            if self.strand == "Forward" and self.paired:
                self.paired = False

            # Write metrics to strand.txt
            if not os.path.exists(self.output):
                os.makedirs(self.output)
            with open(os.path.join(self.output, "strand.txt"), "w") as f:
                f.write("Done setting up the STARsolo run; here are final processing options:\n")
                f.write("=============================================================================\n")
                f.write(f"Sample: {self.snakemake.wildcards.sample}\n")
                f.write(f"Paired-end mode: {self.paired}\n")
                f.write(f"Strand (Forward = 3', Reverse = 5'): {self.strand}, %reads mapped to GeneFull: forward = {PCTFWD} , reverse = {PCTREV}\n")
                f.write(f"CB whitelist: {self.BC}, matches out of 200,000: {nbc3} (v3), {nbc2} (v2), {nbc1} (v1), {nbcA} (multiome)\n")
                f.write(f"CB length: {self.cblen}\n")
                f.write(f"UMI length: {self.umilen}\n")
                f.write(f"GZIP: {self.readcmd}\n")
                f.write("-----------------------------------------------------------------------------\n")
                f.write(f"Read 1 files: {self.fq1}\n")
                f.write("-----------------------------------------------------------------------------\n")
                f.write(f"Read 2 files: {self.fq2}\n")
                f.write("-----------------------------------------------------------------------------\n")

    def run(self):
        with TemporaryDirectory() as tmpdir:
            shell(
                "STAR "
                " --runThreadN {self.snakemake.threads}"
                " --genomeDir {self.index}"
                " {self.input}"
                " {self.readcmd}"
                " {self.out_unmapped}"
                " {self.extra}"
                " --outTmpDir {tmpdir}/STARtmp"
                " --outFileNamePrefix {self.output}/"
                " --soloOutFileNames solo/ features.tsv barcodes.tsv matrix.mtx"
                # " > {self.snakemake.output.aln}"
                " {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)