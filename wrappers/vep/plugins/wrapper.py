# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/10 09:43:54
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import sys
from pathlib import Path
from urllib.request import urlretrieve
from zipfile import ZipFile
from tempfile import NamedTemporaryFile
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        if self.snakemake.log:
            sys.stderr = open(self.snakemake.log[0], "w")

        self.outdir = Path(self.snakemake.output[0])
        self.outdir.mkdir()

    def run(self):
        with NamedTemporaryFile() as tmp:
            urlretrieve(
                "https://github.com/Ensembl/VEP_plugins/archive/release/{release}.zip".format(
                    release=self.snakemake.params.release
                ),
                tmp.name,
            )

            with ZipFile(tmp.name) as f:
                for member in f.infolist():
                    memberpath = Path(member.filename)
                    if len(memberpath.parts) == 1:
                        # skip root dir
                        continue
                    targetpath = self.outdir / memberpath.relative_to(memberpath.parts[0])
                    if member.is_dir():
                        targetpath.mkdir()
                    else:
                        with open(targetpath, "wb") as out:
                            out.write(f.read(member.filename))


if __name__ == '__main__':
    Wrapper(snakemake)