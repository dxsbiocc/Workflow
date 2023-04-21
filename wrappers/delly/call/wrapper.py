__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


from snakemake.shell import shell
from snakemake_wrapper_utils.bcftools import get_bcftools_opts
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.bcftools_opts = get_bcftools_opts(self.snakemake, parse_ref=False, parse_memory=False)

        exclude = self.snakemake.input.get("exclude", "")
        self.exclude = f"-x {exclude}" if exclude else ""
            
    def run(self):
        shell(
            "(OMP_NUM_THREADS={self.snakemake.threads} delly call"
            " -g {self.snakemake.input.ref}"
            " {self.exclude}"
            " {self.extra}"
            " {self.snakemake.input.alns} |"
            # Convert output to specified format
            " bcftools view"
            " {self.bcftools_opts}"
            ") {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)