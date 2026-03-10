# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/09 21:23:43
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file
# ============================================================


import sys
import subprocess
from pathlib import Path
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):

        species = self.snakemake.params.species.lower()
        build = self.snakemake.params.build
        release = int(self.snakemake.params.release)
        out_fmt = Path(self.snakemake.output[0]).suffixes
        self.out_gz = (out_fmt.pop() and True) if out_fmt[-1] == ".gz" else False
        out_fmt = out_fmt.pop().lstrip(".")

        branch = ""
        if release >= 81 and build == "GRCh37":
            # use the special grch37 branch for new releases
            branch = "grch37/"
        elif self.snakemake.params.get("branch"):
            branch = self.snakemake.params.branch + "/"

        flavor = self.snakemake.params.get("flavor", "")
        if flavor:
            flavor += "."

        suffix = ""
        if out_fmt == "gtf":
            suffix = "gtf.gz"
        elif out_fmt == "gff3":
            suffix = "gff3.gz"
        else:
            raise ValueError(
                "invalid format specified. Only 'gtf[.gz]' and 'gff3[.gz]' are currently supported."
            )

        self.url = "ftp://ftp.ensembl.org/pub/{branch}release-{release}/{out_fmt}/{species}/{species_cap}.{build}.{release}.{flavor}{suffix}".format(
            release=release,
            build=build,
            species=species,
            out_fmt=out_fmt,
            species_cap=species.capitalize(),
            suffix=suffix,
            flavor=flavor,
            branch=branch,
        )

    def run(self):
        try:
            if self.out_gz:
                shell(
                    "curl -L {self.url} > {self.snakemake.output[0]} {self.log}")
            else:
                shell(
                    "(curl -L {self.url} | gzip -d > {self.snakemake.output[0]}) {self.log}")
        except subprocess.CalledProcessError as e:
            if self.snakemake.log:
                sys.stderr = open(self.snakemake.log[0], "a")
            print(
                "Unable to download annotation data from Ensembl. "
                "Did you check that this combination of species, build, and release is actually provided?",
                file=sys.stderr,
            )
            exit(1)


if __name__ == '__main__':
    Wrapper(snakemake)
