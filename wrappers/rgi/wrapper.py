import os
import tempfile
from pathlib import Path
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):
    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.input = self.snakemake.input.get("fasta", "")
        self.output = self.snakemake.output

        self.aligner = self.snakemake.params.get("aligner", "DIAMOND")
        self.intype = self.snakemake.params.get("intype", "contig")

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            shell("cut -f 1 -d ' ' {self.input} > {tmpdir}/protein.fa")
            shell(
                "rgi main"
                " -i {tmpdir}/protein.fa"
                " -o {self.output}"
                " -a {self.aligner}"
                " -t {self.intype}"
                " -n {self.snakemake.threads}"
                " {self.extra}"
                " {self.log}"
            )


if __name__ == "__main__":
    wrapper = Wrapper(snakemake)
