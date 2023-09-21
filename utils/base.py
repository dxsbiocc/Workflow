# -*- encoding: utf-8 -*-
# ============================================================
# File        : base.py
# Time        : 2022/08/17 17:15:11
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


import os
from abc import ABCMeta, abstractmethod


class WrapperBase(metaclass=ABCMeta):

    def __init__(self, snakemake) -> None:
        self.snakemake = snakemake
        self.mkdirs()
        # common parameters
        self.extra = snakemake.params.get("extra", "")
        self.tmp = snakemake.params.get("tmp_dir", "/tmp/")
        # log
        self.log = snakemake.log_fmt_shell(stdout=True, stderr=True)
        self.parser()
        self.run()

    def mkdirs(self):
        """Make sure the folder has been created"""
        for file in self.snakemake.output:
            dirname = os.path.dirname(os.path.abspath(file))
            if not os.path.exists(dirname):
                os.makedirs(dirname)
    
    @abstractmethod
    def parser(self):
        """Extract all informations from input、output and parameterss
        """
        pass

    @abstractmethod
    def run(self):
        """Run shell command
        """
        pass