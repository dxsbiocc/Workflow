# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2023/03/10 09:26:18
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import os
from pathlib import Path
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def get_only_child_dir(self, path):
        children = [child for child in path.iterdir() if child.is_dir()]
        assert (
            len(children) == 1
        ), "Invalid VEP cache directory, only a single entry is allowed, make sure that cache was created with the snakemake VEP cache wrapper"
        return children[0]
    
    def parser(self):
        self.fork = "--fork {}".format(self.snakemake.threads) if self.snakemake.threads > 1 else ""
        self.stats = self.snakemake.output.stats
        cache = self.snakemake.input.get("cache", "")
        self.plugins = self.snakemake.input.plugins
        plugin_aux_files = {"LoFtool": "LoFtool_scores.txt", "ExACpLI": "ExACpLI_values.txt"}

        load_plugins = []
        for plugin in self.snakemake.params.plugins:
            if plugin in plugin_aux_files.keys():
                aux_path = os.path.join(self.plugins, plugin_aux_files[plugin])
                load_plugins.append(",".join([plugin, aux_path]))
            else:
                load_plugins.append(",".join([plugin, self.snakemake.input.get(plugin.lower(), "")]))
        self.load_plugins = " ".join(map("--plugin {}".format, load_plugins))

        if self.snakemake.output.calls.endswith(".vcf.gz"):
            self.fmt = "z"
        elif self.snakemake.output.calls.endswith(".bcf"):
            self.fmt = "b"
        else:
            self.fmt = "v"

        fasta = self.snakemake.input.get("fasta", "")
        self.fasta = "--fasta {}".format(fasta) if fasta else ""

        gff = self.snakemake.input.get("gff", "")
        self.gff = "--gff {}".format(gff) if gff else ""

        self.cache = cache
        if cache:
            entrypath = self.get_only_child_dir(self.get_only_child_dir(Path(cache)))
            species = (
                entrypath.parent.name[:-7]
                if entrypath.parent.name.endswith("_refseq")
                else entrypath.parent.name
            )
            release, build = entrypath.name.split("_")
            self.cache = (
                "--offline --cache --dir_cache {cache} --cache_version {release} --species {species} --assembly {build}"
            ).format(cache=cache, release=release, build=build, species=species)

    def run(self):
        shell(
            "(bcftools view '{self.snakemake.input.calls}' | "
            "vep {self.extra} {self.fork} "
            "--format vcf "
            "--vcf "
            "{self.cache} "
            "{self.gff} "
            "{self.fasta} "
            "--dir_plugins {self.plugins} "
            "{self.load_plugins} "
            "--output_file STDOUT "
            "--stats_file {self.stats} | "
            "bcftools view -O{self.fmt} > {self.snakemake.output.calls})"
            " {self.log}"
        )


if __name__ == '__main__':
    Wrapper(snakemake)