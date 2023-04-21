"""Snakemake wrapper for Varscan2 mpileup2indel"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2019, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import os
from snakemake.shell import shell
from snakemake.utils import makedirs
from snakemake_wrapper_utils.java import get_java_opts
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.java_opts = get_java_opts(self.snakemake)
        
        self.pileup = (
            " cat {} ".format(self.snakemake.input[0])
            if not self.snakemake.input[0].endswith("gz")
            else " zcat {} ".format(self.snakemake.input[0])
        )

        # Building output directories
        makedirs(os.path.dirname(self.snakemake.output[0]))

    def run(self):
        shell(
            "varscan mpileup2indel "  # Tool and its subprocess
            "<( {self.pileup} ) "
            "{self.java_opts} {self.extra} "  # Extra parameters
            "> {self.snakemake.output[0]} "  # Path to vcf file
            "{self.log}"  # Logging behaviour
        )


if __name__ == '__main__':
    Wrapper(snakemake)