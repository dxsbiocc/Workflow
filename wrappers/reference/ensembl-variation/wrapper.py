__author__ = "Johannes Köster"
__copyright__ = "Copyright 2019, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import tempfile
import subprocess
import sys
import os
from snakemake.shell import shell
from snakemake.exceptions import WorkflowError
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        species = self.snakemake.params.species.lower()
        release = int(self.snakemake.params.release)
        build = self.snakemake.params.build
        type = self.snakemake.params.type
        chromosome = self.snakemake.params.get("chromosome", "")


        branch = ""
        if release >= 81 and build == "GRCh37":
            # use the special grch37 branch for new releases
            branch = "grch37/"
        elif self.snakemake.params.get("branch"):
            branch = self.snakemake.params.branch + "/"

        if release < 98 and not branch:
            print("Ensembl releases <98 are unsupported.", file=open(self.snakemake.log[0], "w"))
            exit(1)

        if chromosome and type != "all":
            raise ValueError(
                "Parameter chromosome given but chromosome-wise download"
                "is only implemented for type='all'."
            )

        if type == "all":
            if species == "homo_sapiens" and release >= 93:
                chroms = (
                    list(range(1, 23)) + ["X", "Y", "MT"] if not chromosome else [chromosome]
                )
                suffixes = ["-chr{}".format(chrom) for chrom in chroms]
            else:
                if chromosome:
                    raise ValueError(
                        "Parameter chromosome given but chromosome-wise download"
                        "is only implemented for homo_sapiens in releases >=93."
                    )
                suffixes = [""]
        elif type == "somatic":
            suffixes = ["_somatic"]
        elif type == "structural_variations":
            suffixes = ["_structural_variations"]
        else:
            raise ValueError(
                "Unsupported type {} (only all, somatic, structural_variations are allowed)".format(
                    type
                )
            )

        species_filename = species if release >= 91 else species.capitalize()

        self.urls = [
            "ftp://ftp.ensembl.org/pub/{branch}release-{release}/variation/vcf/{species}/{species_filename}{suffix}.{ext}".format(
                release=release,
                species=species,
                suffix=suffix,
                species_filename=species_filename,
                branch=branch,
                ext=ext,
            )
            for suffix in suffixes
            for ext in ["vcf.gz", "vcf.gz.csi"]
        ]
        self.names = [os.path.basename(url) for url in self.urls if url.endswith(".gz")]

    def run(self):
        try:
            gather = "curl {urls}".format(urls=" ".join(map("-O {}".format, self.urls)))
            workdir = os.getcwd()
            with tempfile.TemporaryDirectory() as tmpdir:
                if self.snakemake.input.get("fai"):
                    shell(
                        "(cd {tmpdir}; {gather} && "
                        "bcftools concat -Oz --naive {self.names} > concat.vcf.gz && "
                        "bcftools reheader --fai {workdir}/{self.snakemake.input.fai} concat.vcf.gz "
                        "> {workdir}/{self.snakemake.output}) {self.log}"
                    )
                else:
                    shell(
                        "(cd {tmpdir}; {gather} && "
                        "bcftools concat -Oz --naive {self.names} "
                        "> {workdir}/{self.snakemake.output}) {self.log}"
                    )
        except subprocess.CalledProcessError as e:
            if self.snakemake.log:
                sys.stderr = open(self.snakemake.log[0], "a")
            print(
                "Unable to download variation data from Ensembl. "
                "Did you check that this combination of species, build, and release is actually provided? ",
                file=sys.stderr,
            )
            exit(1)


if __name__ == "__main__":
    Wrapper(snakemake)