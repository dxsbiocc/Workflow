# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/12/07 11:15:46
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):
    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.log = self.snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)
        required_thread_nb = 1

        genome = self.snakemake.input["genome"]
        if genome.endswith(".gz"):
            self.genome = f"<( gzip --stdout --decompress {genome} )"
            required_thread_nb += 1  # Add a thread for gzip uncompression
        elif genome.endswith(".bz2"):
            self.genome = f"<( bzip2 --stdout --decompress {genome} )"
            required_thread_nb += 1  # Add a thread for bzip2 uncompression
        else:
            self.genome = genome

        if self.snakemake.threads < required_thread_nb:
            raise ValueError(
                f"Salmon decoy wrapper requires exactly {required_thread_nb} threads, "
                f"but only {self.snakemake.threads} were provided"
            )

        sequences = [
            self.snakemake.input["transcriptome"],
            self.snakemake.input["genome"],
            self.snakemake.output["gentrome"],
        ]
        if all(fasta.endswith(".gz") for fasta in sequences):
            # Then all input sequences are gzipped. The output will also be gzipped.
            pass
        elif all(fasta.endswith(".bz2") for fasta in sequences):
            # Then all input sequences are bgzipped. The output will also be bgzipped.
            pass
        elif all(fasta.endswith((".fa", ".fna", ".fasta")) for fasta in sequences):
            # Then all input sequences are raw fasta. The output will also be raw fasta.
            pass
        else:
            raise ValueError(
                "Mixed compression status: Either all fasta sequences are compressed "
                "with the *same* compression algorithm, or none of them are compressed."
            )
    
    def run(self):
        # Gathering decoy sequences names
        # Sed command works as follow:
        # -n       = do not print all lines
        # s/ .*//g = Remove anything after spaces. (remove comments)
        # s/>//p  = Remove '>' character at the begining of sequence names. Print names.
        shell("( sed -n 's/ .*//g;s/>//p' {self.genome} ) > {self.snakemake.output.decoys} {self.log}")

        # Building big gentrome file
        shell(
            "cat {self.snakemake.input.transcriptome} {self.snakemake.input.genome} "
            "> {self.snakemake.output.gentrome} {self.log}"
        )
    

if __name__ == '__main__':
    Wrapper(snakemake)
