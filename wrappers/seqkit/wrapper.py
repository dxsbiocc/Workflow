from snakemake.shell import shell
from snakemake_wrapper_utils.base import WrapperBase


class Wrapper(WrapperBase):
    def __init__(self, snakemake) -> None:
        super().__init__(snakemake)

    def parser(self):
        commond = self.snakemake.params.get("command", "seq")
        assert commond in [
            "amplicon", "bam", "common", "concat", "convert", "duplicate",
            "fa2fq", "faidx", "fish", "fq2fa", "fx2tab", "genautocomplete",
            "grep", "head", "head-genome", "locate", "merge-slides", "mutate",
            "pair", "range", "rename", "replace", "restart", "rmdup",
            "sample", "sana", "scat", "seq", "shuffle", "sliding", "sort", 
            "split", "split2", "stats", "subseq", "sum", "tab2fx", "translate"
        ], "command not support!"

        self.extra_input = " ".join(
            [
                f"--{key.replace('_','-')} {value}"
                if key in ["bed", "gtf"]
                else f"--{key.replace('_','-')}-file {value}"
                for key, value in self.snakemake.input.items()
            ][1:]
        )

        self.extra_output = " ".join(
            [
                f"--{key.replace('_','-')} {value}"
                if key in ["read1", "read2"]
                else f"--{key.replace('_','-')}-file {value}"
                for key, value in self.snakemake.output.items()
            ][1:]
        )

    
    def run(self):
        shell(
            "seqkit {self.snakemake.params.command}"
            " --threads {self.snakemake.threads}"
            " {self.extra_input}"
            " {self.extra_output}"
            " {self.extra}"
            " --out-file {self.snakemake.output[0]}"
            " {self.snakemake.input[0]}"
            " {self.log}"
        )
    

if __name__ == "__main__":
    wrapper = Wrapper(snakemake)
