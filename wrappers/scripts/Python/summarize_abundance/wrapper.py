# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2024/01/11 20:01:59
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : Summarize abundance based on group mapping
# ============================================================


import sys
from pathlib import Path
import pandas as pd
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        # Input files
        self.input_file = self.snakemake.input.get("abundance")
        self.map_file = self.snakemake.input.get("map")

        if not self.input_file:
            raise ValueError(
                "Input abundance file must be specified in snakemake.input.abundance"
            )
        if not self.map_file:
            raise ValueError("Input map file must be specified in snakemake.input.map")

        # Output prefix - use params if specified, otherwise extract from output directory
        # The script generates multiple output files: {prefix}.{grp}.{exprType}.txt
        if self.snakemake.params.get("output_prefix"):
            self.output_prefix = self.snakemake.params.get("output_prefix")
        else:
            # Try to infer from output files
            # Get output directory from first output file
            output_files = self.snakemake.output
            if isinstance(output_files, dict):
                first_output = list(output_files.values())[0]
            else:
                first_output = output_files

            output_path = Path(first_output)
            # Use directory + a default prefix name
            # Or extract from filename if it follows the pattern
            stem = output_path.stem
            parts = stem.split(".")
            if len(parts) >= 3:
                # Format: prefix.{grp}.{exprType}
                base_name = ".".join(parts[:-2])
                self.output_prefix = str(output_path.parent / base_name)
            else:
                # Use directory + default name
                self.output_prefix = str(output_path.parent / "summarized")

        # Parameters with defaults matching the script
        self.grp_column = self.snakemake.params.get("grp_column", "2")
        self.grp_sep = self.snakemake.params.get("grp_sep", ",")
        self.subitem_keep = self.snakemake.params.get("subitem_keep", "all")
        self.abundance_keep = self.snakemake.params.get("abundance_keep", "sum")
        self.norm_type = self.snakemake.params.get("norm_type", "raw")
        self.dropkeycol = self.snakemake.params.get("dropkeycol", False)
        self.verbose = self.snakemake.params.get("verbose", False)
        self.debug = self.snakemake.params.get("debug", False)

    def norm_expr(self, probe_expr_df, norm_type):
        """Normalize expression data"""
        if norm_type == "raw":
            return probe_expr_df
        elif norm_type == "cpm":
            norm_factor = probe_expr_df.sum()
            return probe_expr_df * 1000000 / norm_factor
        else:
            raise ValueError(f"Unsupported normalization type: {norm_type}")

    def run(self):
        # Parse parameters
        grp_columnL = [int(i) - 1 for i in self.grp_column.split(",")]
        grp_column_len = len(grp_columnL)
        grp_sepL = [i for i in self.grp_sep.split("+")]
        grp_sep_len = len(grp_sepL)
        norm_typeL = [i for i in self.norm_type.split(",")]

        # Handle separator list
        if grp_sep_len == 1:
            grp_sepL = grp_sepL * grp_column_len
        else:
            if grp_sep_len != grp_column_len:
                raise ValueError("Unequal number of columns (-c) and separators (-s)")

        if self.debug:
            print(f"grp_columnL={grp_columnL}", file=sys.stderr)
            print(f"grp_sepL={grp_sepL}", file=sys.stderr)

        # Read abundance matrix
        probe_expr_df = pd.read_csv(self.input_file, sep="\t", header=0, index_col=0)
        probe_expr_df.index = probe_expr_df.index.map(str)

        if self.debug:
            print(
                f"Probes example in expression file: {list(probe_expr_df.index[:10])}.",
                file=sys.stderr,
            )
            print(
                f"Read in probe expression matrix: {probe_expr_df.shape}.",
                file=sys.stderr,
            )

        # Normalize expression data
        abundanceL = []
        for norm_type in norm_typeL:
            abundanceL.append(self.norm_expr(probe_expr_df, norm_type))

        # Read map file
        usecols = [0]
        usecols.extend(grp_columnL)
        probe_map_df = pd.read_csv(
            self.map_file,
            sep="\t",
            header=0,
            index_col=None,
            dtype=str,
            usecols=usecols,
        )

        if self.debug:
            print(f"Read in probe-map matrix: {probe_map_df.shape}.", file=sys.stderr)
            print(f"Read in probe-map matrix:\n{probe_map_df.head()}.", file=sys.stderr)

        keyName = list(probe_map_df.columns)[0]
        grp_column_names = list(probe_map_df.columns)[1:]

        # Process each group column
        for grp, sep in zip(grp_column_names, grp_sepL):
            # Filter out empty group values
            grpDF = probe_map_df.loc[probe_map_df[grp] != "", [keyName, grp]]

            # Process each normalization type
            for exprDF, exprType in zip(abundanceL, norm_typeL):
                # Merge abundance data with group mapping
                grpMergeDF = grpDF.merge(
                    exprDF, left_on=keyName, right_index=True, how="inner"
                )

                # Drop key column if requested
                if self.dropkeycol:
                    grpMergeDF.drop(keyName, axis=1, inplace=True)

                if self.debug:
                    print(f"Merged matrix: {grpMergeDF.shape}.", file=sys.stderr)
                    print(f"Merged matrix:\n{grpMergeDF.head()}.", file=sys.stderr)

                # Split group values by separator
                if sep == "*":
                    sep = ""
                grpMergeDF[grp] = grpMergeDF[grp].str.split(sep)
                grpMergeDF_explode = grpMergeDF.explode(grp)

                if self.debug:
                    print(
                        f"Exploded matrix:\n{grpMergeDF_explode.head()}.",
                        file=sys.stderr,
                    )

                # Filter out empty group values after explode
                grpMergeDF_explode = grpMergeDF_explode.loc[
                    grpMergeDF_explode[grp] != "",
                ]

                # Aggregate by group
                final_abundanceDF = grpMergeDF_explode.groupby(grp, as_index=False).agg(
                    self.abundance_keep
                )

                # Write output file
                filename = f"{self.output_prefix}.{grp}.{exprType}.txt"
                final_abundanceDF.to_csv(filename, sep="\t", index=False)

                if self.verbose:
                    print(f"Output written to: {filename}", file=sys.stderr)


if __name__ == "__main__":
    Wrapper(snakemake)
