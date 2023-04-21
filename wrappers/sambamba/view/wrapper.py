__author__ = "Jan Forster"
__copyright__ = "Copyright 2021, Jan Forster"
__email__ = "j.forster@dkfz.de"
__license__ = "MIT"


import os

from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        in_file = self.snakemake.input[0]
        if in_file.endswith(".sam") and ("-S" not in self.extra or "--sam-input" not in self.extra):
            self.extra += " --sam-input"

    def run(self):
        shell(
            "sambamba view"
            " {self.extra}"
            " -t {self.snakemake.threads}"
            " {self.snakemake.input[0]}"
            " > {self.snakemake.output[0]}"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)