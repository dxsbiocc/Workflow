__author__ = "Johannes Köster"
__copyright__ = "Copyright 2019, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import subprocess as sp
import sys
from itertools import product
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        species = self.snakemake.params.species.lower()
        release = int(self.snakemake.params.release)
        build = self.snakemake.params.build

        branch = ""
        if release >= 81 and build == "GRCh37":
            # use the special grch37 branch for new releases
            branch = "grch37/"
        elif self.snakemake.params.get("branch"):
            branch = self.snakemake.params.branch + "/"

        spec = ("{build}" if int(release) > 75 else "{build}.{release}").format(
            build=build, release=release
        )

        self.suffixes = ""
        datatype = self.snakemake.params.get("datatype", "")
        chromosome = self.snakemake.params.get("chromosome", "")
        if datatype == "dna":
            if chromosome:
                self.suffixes = ["dna.chromosome.{}.fa.gz".format(chromosome)]
            else:
                self.suffixes = ["dna.primary_assembly.fa.gz", "dna.toplevel.fa.gz"]
        elif datatype == "cdna":
            self.suffixes = ["cdna.all.fa.gz"]
        elif datatype == "cds":
            self.suffixes = ["cds.all.fa.gz"]
        elif datatype == "ncrna":
            self.suffixes = ["ncrna.fa.gz"]
        elif datatype == "pep":
            self.suffixes = ["pep.all.fa.gz"]
        else:
            raise ValueError(
                "invalid datatype, must be one of dna, cdna, cds, ncrna, pep")

        if chromosome:
            if not datatype == "dna":
                raise ValueError(
                    "invalid datatype, to select a single chromosome the datatype must be dna"
                )

        spec = spec.format(build=build, release=release)
        self.url_prefix = f"ftp://ftp.ensembl.org/pub/{branch}release-{release}/fasta/{species}/{datatype}/{species.capitalize()}.{spec}"

    def run(self):
        success = False
        for suffix in self.suffixes:
            url = f"{self.url_prefix}.{suffix}"

            try:
                shell("curl -sSf {url} > /dev/null 2> /dev/null")
            except sp.CalledProcessError:
                continue

            shell("(curl -L {url} | gzip -d > {self.snakemake.output[0]}) {self.log}")
            success = True
            break

        if not success:
            if len(self.suffixes) > 1:
                url = f"{self.url_prefix}.[{'|'.join(self.suffixes)}]"
            else:
                url = f"{self.url_prefix}.{self.suffixes[0]}"
            print(
                f"Unable to download requested sequence data from Ensembl ({url}). "
                "Please check whether above URL is currently available (might be a temporal server issue). "
                "Apart from that, did you check that this combination of species, build, and release is actually provided?",
                file=sys.stderr,
            )
            exit(1)


if __name__ == '__main__':
    Wrapper(snakemake)