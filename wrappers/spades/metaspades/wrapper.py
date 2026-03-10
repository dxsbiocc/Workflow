"""Snakemake wrapper for metaspades."""

import os
import shutil
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):
    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        # infer output directory
        if hasattr(self.snakemake.output, "dir"):
            self.output_dir = self.snakemake.output.dir

        else:
            # get output_dir file from output
            if hasattr(self.snakemake.output, "contigs"):
                output_file = self.snakemake.output.contigs
            elif hasattr(self.snakemake.output, "scaffolds"):
                output_file = self.snakemake.output.scaffolds
            else:
                output_file = self.snakemake.output[0]

            self.output_dir = os.path.split(output_file)[0]

        # parse params
        self.kmers = self.snakemake.params.get("k", "'auto'")
        self.mode = self.snakemake.params.get("mode", "meta")

        assert self.mode in [
            "meta",
            "isolate",
            "sc",
            "plasmid",
            "metaplasmid",
            "metaviral",
            "rna",
            "rnaviral",
            "bio",
            "corona",
            "sewage",
        ], "Invalid mode: {}".format(self.mode)

        if hasattr(self.snakemake.resources, "mem_mb"):
            mem_gb = self.snakemake.resources.mem_mb // 1000
            self.memory_requirements = f" --memory {mem_gb}"
        else:
            self.memory_requirements = ""

        if not os.path.exists(os.path.join(self.output_dir, "params.txt")):
            # parse short reads
            if hasattr(self.snakemake.input, "reads"):
                reads = self.snakemake.input.reads

                assert (
                    len(reads) > 1
                ), "Metaspades needs a paired end library. This means you should supply at least 2 fastq files in the rule input."

                assert (
                    type(reads[0]) == str
                ), f"Metaspades allows only 1 library. Therefore reads need to be strings got {reads}"

                input_arg = " --pe1-1 {0} --pe1-2 {1} ".format(*reads)

                if len(reads) >= 3:
                    input_arg += " --pe1-m {2}".format(*reads)

                    if len(reads) >= 4:
                        input_arg += " --pe1-s {3}".format(*reads)
            elif hasattr(self.snakemake.input, "r1") and hasattr(
                self.snakemake.input, "r2"
            ):
                if isinstance(self.snakemake.input.r1, list):
                    r1 = ",".join(self.snakemake.input.r1)
                else:
                    r1 = self.snakemake.input.r1
                if isinstance(self.snakemake.input.r2, list):
                    r2 = ",".join(self.snakemake.input.r2)
                else:
                    r2 = self.snakemake.input.r2
                # check if the number of R1 and R2 files are the same
                if len(self.snakemake.input.r1) != len(self.snakemake.input.r2):
                    raise ValueError("The number of R1 and R2 files must be the same.")
                input_arg = f" -1 {r1} -2 {r2} "
            elif hasattr(self.snakemake.input, "unpaired"):
                input_arg = f" -s {self.snakemake.input.unpaired} "

            # parse long reads
            for longread_name in ["pacbio", "nanopore"]:
                if hasattr(self.snakemake.input, longread_name):
                    input_arg += " --{name} {}".format(
                        name=longread_name, **self.snakemake.input
                    )
            self.input_arg = input_arg

    def run(self):
        if not os.path.exists(os.path.join(self.output_dir, "params.txt")):
            shell(
                "spades.py"
                " --{self.mode} "
                " --threads {self.snakemake.threads} "
                " {self.memory_requirements} "
                " -o {self.output_dir} "
                " -k {self.kmers} "
                " {self.input_arg} "
                " {self.extra} "
                " > {self.log[0]} 2>&1 "
            )

        else:
            # params.txt file exitst already I restart from previous run

            shell(
                "echo '\n\nRestart Spades \n Remove pipline_state file copy files to force copy files if necessary.' >> {log[0]}"
            )

            shell(
                "rm -f {self.output_dir}/pipeline_state/stage_*_copy_files 2>> {self.log}"
            )

            shell(
                "spades.py --meta "
                " --restart-from last "
                " --threads {self.snakemake.threads} "
                " {self.memory_requirements} "
                " -o {self.output_dir} "
                " >> {self.log[0]} 2>&1 "
            )

        # Rename/ move output files
        Output_key_mapping = {
            "contigs": "contigs.fasta",
            "scaffolds": "scaffolds.fasta",
            "graph": "assembly_graph_with_scaffolds.gfa",
        }

        has_named_output = False
        for key in Output_key_mapping:
            if hasattr(self.snakemake.output, key):
                has_named_output = True
                file_produced = os.path.join(self.output_dir, Output_key_mapping[key])
                file_renamed = getattr(self.snakemake.output, key)

                if file_produced != file_renamed:
                    shutil.move(file_produced, file_renamed)

        if not has_named_output:
            file_produced = os.path.join(self.output_dir, "contigs.fasta")
            file_renamed = self.snakemake.output[0]

            if file_produced != file_renamed:
                shutil.move(file_produced, file_renamed)


if __name__ == "__main__":
    Wrapper(snakemake)
