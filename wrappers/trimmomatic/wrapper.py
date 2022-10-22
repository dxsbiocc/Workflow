# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2022/10/20 10:43:35
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        self.java_opts = get_java_opts(snakemake)
        super().__init__(snakemake)
    
    def parser(self):
        compression_level = self.snakemake.params.get("compression_level", "-5")
        self.trimmer = " ".join(self.snakemake.params.trimmer)

        # Distribute threads
        num = len(self.snakemake.input)
        if num == 1:
            input_files, output_files = [self.snakemake.input[0]], [self.snakemake.output[0]]
            self.type = 'SE'
        else:
            input_files = [self.snakemake.input.fq1, self.snakemake.input.fq2]
            output_files = [
                self.snakemake.output.fq1,
                self.snakemake.output.fq1_unpaired,
                self.snakemake.output.fq2,
                self.snakemake.output.fq2_unpaired,
            ]
            self.type = 'PE'
        self.inputs = input_files
        self.outputs = output_files
        # self.trimmomatic_threads, input_threads, output_threads = self.distribute_threads(
        #     input_files, output_files, self.snakemake.threads
        # )

        # self.inputs = " ".join([
        #     self.compose_input_gz(filename, input_threads) for filename in input_files
        # ])

        # self.outputs = " ".join([
        #     self.compose_output_gz(filename, output_threads, compression_level)
        #     for filename in output_files
        # ])

    def run(self):
        shell(
            "trimmomatic {self.type} -threads {self.trimmomatic_threads} {self.java_opts} {self.extra} "
            "{self.inputs} {self.outputs} "
            "{self.trimmer} "
            "{self.log}"
        )

    # Distribute available threads between trimmomatic itself and any potential pigz instances
    def distribute_threads(self, input_files, output_files, available_threads):
        gzipped_input_files = sum(1 for file in input_files if file.endswith(".gz"))
        gzipped_output_files = sum(1 for file in output_files if file.endswith(".gz"))
        potential_threads_per_process = available_threads // (
            1 + gzipped_input_files + gzipped_output_files
        )
        if potential_threads_per_process > 0:
            # decompressing pigz creates at most 4 threads
            pigz_input_threads = (
                min(4, potential_threads_per_process) if gzipped_input_files != 0 else 0
            )
            pigz_output_threads = (
                (available_threads - pigz_input_threads * gzipped_input_files)
                // (1 + gzipped_output_files)
                if gzipped_output_files != 0
                else 0
            )
            trimmomatic_threads = (
                available_threads
                - pigz_input_threads * gzipped_input_files
                - pigz_output_threads * gzipped_output_files
            )
        else:
            # not enough threads for pigz
            pigz_input_threads = 0
            pigz_output_threads = 0
            trimmomatic_threads = available_threads
        return trimmomatic_threads, pigz_input_threads, pigz_output_threads


    def compose_input_gz(self, filename, threads):
        if filename.endswith(".gz") and threads > 0:
            return "<(pigz -p {threads} --decompress --stdout {filename})".format(
                threads=threads, filename=filename
            )
        return filename


    def compose_output_gz(self, filename, threads, compression_level):
        if filename.endswith(".gz") and threads > 0:
            return ">(pigz -p {threads} {compression_level} > {filename})".format(
                threads=threads, compression_level=compression_level, filename=filename
            )
        return filename


if __name__ == '__main__':
    Wrapper(snakemake)