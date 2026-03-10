import os
import tempfile
from pathlib import Path
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):
    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.input_file = self.snakemake.input.get("input_file", "")
        self.output_file = self.snakemake.output.get("output_file", "")

        trans_file = self.snakemake.output.get("trans_file", "")
        if trans_file:
            self.extra += f" -a {trans_file}"
        nuc_file = self.snakemake.output.get("nuc_file", "")
        if nuc_file:
            self.extra += f" -d {nuc_file}"
        start_file = self.snakemake.output.get("start_file", "")
        if start_file:
            self.extra += f" -s {start_file}"

        training_file = self.snakemake.input.get("training_file", "")
        if not os.path.exists(training_file):
            training_file = self.snakemake.output.get("training_file", "")
            if training_file:
                self.extra += f" -g {training_file}"

    def run(self):
        shell(
            "prodigal"
            " -i {self.input_file}"
            " -o {self.output_file}"
            " {self.extra}"
            " {self.log}"
        )


if __name__ == "__main__":
    wrapper = Wrapper(snakemake)
