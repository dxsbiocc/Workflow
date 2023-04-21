__author__ = "Johannes Köster, Patrik Smeds"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import os
import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        self.expected_endings = [
            "INT",
            "D",
            "SI",
            "INV",
            "INV_final",
            "TD",
            "LI",
            "BP",
            "CloseEndMapped",
            "RP",
        ]


    def split_file_name(self, file_parts, file_ending_index):
        return (
            "_".join(file_parts[:file_ending_index]),
            "_".join(file_parts[file_ending_index:]),
        )


    def process_input_path(self, input_file):
        """
        :params input_file: Input file from rule, ex /path/to/file/all_D or /path/to/file/all_INV_final
        :return: ""/path/to/file", "all"

        """
        file_path, file_name = os.path.split(input_file)
        file_parts = file_name.split("_")
        # seperate ending and name, to name: all ending: D or name: all ending: INV_final
        file_name, file_ending = self.split_file_name(
            file_parts, -2 if file_name.endswith("_final") else -1
        )
        if not file_ending in self.expected_endings:
            raise Exception("Unexpected variant type: " + file_ending)
        return file_path, file_name

    def run(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            self.input_flag = "-p"
            self.input_file = self.snakemake.input.get("pindel")
            if isinstance(self.input_file, list) and len(self.input_file) > 1:
                self.input_flag = "-P"
                input_path, input_name = self.process_input_path(self.input_file[0])
                self.input_file = os.path.join(input_path, input_name)
                for variant_input in self.snakemake.input.pindel:
                    if not variant_input.startswith(self.input_file):
                        raise Exception(
                            "Unable to extract common path from multi file input, expect path is: "
                            + self.input_file
                        )
                    if not os.path.isfile(variant_input):
                        raise Exception('Input "' + self.input_file + '" is not a file!')
                    os.symlink(
                        os.path.abspath(variant_input),
                        os.path.join(tmpdirname, os.path.basename(variant_input)),
                    )
                self.input_file = os.path.join(tmpdirname, input_name)

            shell(
                "pindel2vcf"
                " {self.extra}"
                " {self.input_flag}"
                " {self.input_file}"
                " -r {self.snakemake.input.ref}"
                " -R {self.snakemake.params.refname}"
                " -d {self.snakemake.params.refdate}"
                " -v {self.snakemake.output[0]}"
                " {self.log}"
            )


if __name__ == '__main__':
    Wrapper(snakemake)