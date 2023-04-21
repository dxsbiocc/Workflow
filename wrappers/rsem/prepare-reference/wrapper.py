__author__ = "Brett Copeland"
__copyright__ = "Copyright 2021, Brett Copeland"
__email__ = "brcopeland@ucsd.edu"
__license__ = "MIT"


import os

from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        # the reference_name argument is inferred by stripping the .seq suffix from
        # the output.seq value
        output_directory = os.path.dirname(os.path.abspath(self.snakemake.output.seq))
        seq_file = os.path.basename(self.snakemake.output.seq)
        if seq_file.endswith(".seq"):
            self.reference_name = os.path.join(output_directory, seq_file[:-4])
        else:
            raise Exception("output.seq has an invalid file suffix (must be .seq)")

        for output_variable, output_path in self.snakemake.output.items():
            if not os.path.abspath(output_path).startswith(self.reference_name):
                raise Exception(
                    "the path for {} is inconsistent with that of output.seq".format(output_variable)
                )

    def run(self):
        shell(
            "rsem-prepare-reference"
            " --num-threads {self.snakemake.threads}"
            " {self.extra}"
            " {self.snakemake.input.reference_genome}"
            " {self.reference_name} "
            "{self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)