import os
import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):
    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        
        ref = self.snakemake.input.get("ref", "")
        if ref:
            self.extra += f" -r {ref}"

        gff = self.snakemake.input.get("gff", "")
        if gff:
            self.extra += f" --features {gff}"

        pe1 = self.snakemake.input.get("pe1", "")
        if pe1:
            self.extra += f" --pe1 {pe1}"
        pe2 = self.snakemake.input.get("pe2", "")
        if pe2:
            self.extra += f" --pe2 {pe2}"
        pe12 = self.snakemake.input.get("pe12", "")
        if pe12:
            self.extra += f" --pe12 {pe12}"
        mp1 = self.snakemake.input.get("mp1", "")
        if mp1:
            self.extra += f" --mp1 {mp1}"
        mp2 = self.snakemake.input.get("mp2", "")
        if mp2:
            self.extra += f" --mp2 {mp2}"
        mp12 = self.snakemake.input.get("mp12", "")
        if mp12:
            self.extra += f" --mp12 {mp12}"
        single = self.snakemake.input.get("single", "")
        if single:
            self.extra += f" --single {single}"
        pacbio = self.snakemake.input.get("pacbio", "")
        if pacbio:
            self.extra += f" --pacbio {pacbio}"
        nanopore = self.snakemake.input.get("nanopore", "")
        if nanopore:
            self.extra += f" --nanopore {nanopore}"
        ref_bam = self.snakemake.input.get("ref_bam", "")
        if ref_bam:
            self.extra += f" --ref-bam {ref_bam}"
        ref_sam = self.snakemake.input.get("ref_sam", "")
        if ref_sam:
            self.extra += f" --ref-sam {ref_sam}"
        bam = self.snakemake.input.get("bam", "")
        if bam:
            if isinstance(bam, list):
                bam = ",".join(bam)
            self.extra += f" --bam {bam}"
        sam = self.snakemake.input.get("sam", "")
        if sam:
            if isinstance(sam, list):
                sam = ",".join(sam)
            self.extra += f" --sam {sam}"
        sv_bedpe = self.snakemake.input.get("sv_bedpe", "")
        if sv_bedpe:
            self.extra += f" --sv-bedpe {sv_bedpe}"

    
    def run(self):
        shell(
            "(quast --threads {snakemake.threads}"
            " {self.extra}"
            " -o {self.snakemake.output}"
            " {snakemake.input.fasta}"
            ") {self.log}"
        )


if __name__ == "__main__":
    wrapper = Wrapper(snakemake)
