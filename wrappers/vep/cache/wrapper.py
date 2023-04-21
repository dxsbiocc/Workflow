# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/10 09:43:37
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        try:
            self.release = int(self.snakemake.params.release)
        except ValueError:
            raise ValueError("The parameter release is supposed to be an integer.")

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            # We download the cache tarball manually because vep_install does not consider proxy settings (in contrast to curl).
            # See https://github.com/bcbio/bcbio-nextgen/issues/1080
            vep_dir = "vep" if self.release >= 97 else "VEP"
            cache_tarball = (
                f"{self.snakemake.params.species}_vep_{self.release}_{self.snakemake.params.build}.tar.gz"
            )
            log = self.snakemake.log_fmt_shell(stdout=True, stderr=True)
            shell(
                "curl -L ftp://ftp.ensembl.org/pub/release-{self.snakemake.params.release}/"
                "variation/{vep_dir}/{cache_tarball} "
                "-o {tmpdir}/{cache_tarball} {log}"
            )

            log = self.snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
            shell(
                "vep_install --AUTO c "
                "--SPECIES {self.snakemake.params.species} "
                "--ASSEMBLY {self.snakemake.params.build} "
                "--VERSION {self.release} "
                "--CACHEURL {tmpdir} "
                "--CACHEDIR {self.snakemake.output} "
                "--CONVERT "
                "--NO_UPDATE "
                "{self.extra} {log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)