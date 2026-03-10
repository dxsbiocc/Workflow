# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2024/01/11 20:01:59
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file
# ============================================================


import os
from snakemake import shell
from snakemake_wrapper_utils.base import WrapperBase
from snakemake_wrapper_utils.snakemake import get_mem


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        self.checkm_db_url = "https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz"
        super().__init__(snakemake)

    def parser(self):
        self.bins_dir = self.snakemake.input.get("bins_dir")
        bins = self.snakemake.params.get("bins")
        if isinstance(bins, str):
            bins = [bins]

        LETTERS = "ABC"
        for i, bin in enumerate(bins[:3]):  # only take the first 3 bins
            self.extra += f" -{LETTERS[i]} {self.bins_dir}/{bin}_bins"

        self.output = self.snakemake.output
        self.memory_requirements = get_mem(self.snakemake, out_unit="MiB")

        checkm_db = str(self.snakemake.params.get("checkm_db", ""))
        if not checkm_db and "skip-checkm" not in self.extra:
            self.extra += " --skip-checkm"

        if not "skip-checkm" in self.extra:
            if not os.path.exists(checkm_db):
                os.makedirs(checkm_db)
                shell(
                    f"wget -O {os.path.join(checkm_db, 'checkm.tar.gz')} {self.checkm_db_url}"
                )
                shell(
                    f"tar -xvf {os.path.join(checkm_db, 'checkm.tar.gz')} -C {checkm_db}"
                )
            # set the checkm database path
            try:
                shell(f"printf '{checkm_db}' | checkm data setRoot")
            except:
                try:
                    print(
                        f"CheckM database path {checkm_db} setup failed! Trying again..."
                    )
                    shell(f"printf '{checkm_db}' | checkm data setRoot")
                except:
                    raise ValueError(
                        f"CheckM database path {checkm_db} setup failed! Please check if checkm is installed correctly."
                    )
            else:
                print(f"CheckM database path {checkm_db} setup successfully!")

    def run(self):
        shell(
            "metawrap bin_refinement"
            " -o {self.output}"
            " -t {self.snakemake.threads}"
            " -m {self.memory_requirements}"
            " {self.extra}"
            " {self.log}"
        )


if __name__ == "__main__":
    Wrapper(snakemake)
