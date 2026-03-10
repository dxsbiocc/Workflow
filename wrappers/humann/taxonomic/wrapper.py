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
import sys
from tempfile import TemporaryDirectory
from snakemake import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.input_dir = os.path.commonpath(self.snakemake.input)
        # output
        self.output_tsv = self.snakemake.output.tsv
        self.output_spf = self.snakemake.output.spf

    def run(self):
        # merge tsv
        shell(
            "merge_metaphlan_tables.py"
            " {self.input_dir}/*/*_metaphlan_profile.tsv"
            " | sed 's/_1_metaphlan//g'"
            " | tail -n+2"
            " | sed '1 s/clade_name/ID/'"
            " > {self.output_tsv}"
        )

        # generate spf
        temp_spf = self.output_spf + ".tmp"
        self.metaphlan_to_stamp(self.output_tsv, temp_spf)
        # sort and remove duplicates (matching original Perl script behavior)
        shell("sort -r {temp_spf} | uniq > {self.output_spf} && rm {temp_spf}")
        # filter unclassified
        shell(
            "grep -v 'unclassified'"
            " {self.output_spf}"
            f" > {self.output_spf.replace('.spf', '.classified.spf')}"
        )

    def metaphlan_to_stamp(self, input_file, output_file):
        """
        Convert MetaPhlAn output to STAMP format.
        This method replicates the functionality of metaphlan_to_stamp.pl
        """
        # Default ranks
        taxa_ranks = [
            "Kingdom",
            "Phylum",
            "Class",
            "Order",
            "Family",
            "Genus",
            "Species",
        ]

        try:
            with open(input_file, "r") as f_in, open(output_file, "w") as f_out:
                # Read header
                header = f_in.readline().strip()
                cols = header.split("\t")

                # Keep sample names (skip the first 'ID' column)
                samples = cols[1:]

                # Check for MetaPhlAn 2 indicator
                # Start with assumption that this is not MetaPhlAn2
                metaphlan2_version = False
                printed_header = False

                # Process lines
                for line in f_in:
                    line = line.strip()

                    # Check if line starts with '#' (MetaPhlAn2 indicator)
                    if line.startswith("#"):
                        metaphlan2_version = True
                        # MetaPhlAn2 goes down to the strain level
                        if "Strain" not in taxa_ranks:
                            taxa_ranks.append("Strain")
                        continue

                    # Print header once
                    if not printed_header:
                        f_out.write("\t".join(taxa_ranks + samples) + "\n")
                        printed_header = True

                    # Split data
                    data = line.split("\t")
                    if not data:
                        continue

                    # Split taxonomy (first column) by pipe
                    taxonomy_str = data[0]
                    abundance_data = data[1:]

                    taxonomy = taxonomy_str.split("|")

                    unclassified_flag = False
                    final_taxonomy = []
                    should_skip = False

                    # Process each rank
                    for x in range(len(taxa_ranks)):
                        if x < len(taxonomy):
                            # The level exists in the data
                            level_val = taxonomy[x]
                            if "unclassified" in level_val:
                                level_val = "unclassified"
                                unclassified_flag = True
                            final_taxonomy.append(level_val)
                        elif unclassified_flag:
                            # Parent was unclassified, fill down
                            final_taxonomy.append("unclassified")
                        else:
                            # Not defined and not an unclassified parent
                            # Do not include in output
                            should_skip = True
                            break

                    if not should_skip:
                        f_out.write("\t".join(final_taxonomy + abundance_data) + "\n")

                # If header was never printed (empty file), print it now
                if not printed_header:
                    f_out.write("\t".join(taxa_ranks + samples) + "\n")

        except FileNotFoundError:
            print(f"Error: File '{input_file}' not found.", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"Error processing file: {e}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    Wrapper(snakemake)
