# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2024/01/06 11:20:56
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import os, sys
import pandas as pd
from glob import glob
from snakemake_wrapper_utils.base import WrapperBase, get_logger


class Wrapper(WrapperBase):
    
    def __init__(self, snakemake) -> None:
        self.logger = get_logger(snakemake.rule, filename=str(snakemake.log))
        super().__init__(snakemake)

    def parser(self):
        # quant file path
        self.input = self.snakemake.input.get("quant")
        assert self.input, "The input quant file is not specified!"
        # output file
        commonpath = os.path.commonpath(self.input)
        self.output = self.snakemake.output.get('exp', os.path.join(commonpath, "expression.tsv"))

        quant = self.snakemake.params.get("quantifier")
        # load annotation
        anno_file = self.snakemake.input.get("gtf")
        self.logger.info("parse GTF file ...")
        try:
            self.gtf_parser(anno_file, quant)
        except Exception as e:
            self.logger.error(e)
            sys.exit(1)
        
        self.logger.info(f"merge {quant} result ...")
        try:
            if quant == "salmon":
                self.salmon_parser()
            elif quant == "kallisto":
                self.kallisto_parser()
            elif quant == "rsem":
                self.rsem_parser()
            elif quant == "featurecounts":
                self.featurecounts_parser()
            elif quant == "htseq":
                self.htseq_parser()
            else:
                raise ValueError(f"The quantifier `{quant}` is not supported!")
        except Exception as e:
            self.logger.error(e)
            sys.exit(1)
    
    def helper(self, group):
        """Helper function to compute gene length (all exon, exclude intron)
        """
        merged = []
        for _, row in group.sort_values('start').iterrows():
            if not merged or merged[-1][1] < row['start']:
                merged.append([row['start'], row['end']])
            else:
                merged[-1][1] = max(merged[-1][1], row['end'])
        return sum([i[1] - i[0] for i in merged])
    
    def gtf_parser(self, anno_file, quant):
        gtf = pd.read_csv(anno_file, skiprows=5, sep="\t", usecols=[0, 2, 3, 4, 8], 
                            header=None, names=['chrom', 'type', 'start', 'end', 'attr'])
        if quant in ["rsem", "htseq", "featurecounts"]:
            attr = gtf.query('type == "gene"')['attr']  # gencode
            if attr.empty:
                attr = gtf.query('type == "transcript"')['attr']  # refGene
            info = attr.str.extractall(
                r'gene_id "(?P<gene_id>.*?)";.*gene_name "(?P<gene_name>.*?)";').set_index('gene_id')
            # dupliactes gene name
            dup_genes = info.loc[info.duplicated(), 'gene_name'].unique()
            # remove duplicates
            self.id2gene = info.query('gene_name not in @dup_genes')
        elif quant in ["salmon", "kallisto"]:
            attr = gtf.query('type == "transcript"')['attr']
            self.trans2gene = attr.str.extractall(
                r'transcript_id "(?P<transcript_id>.*?)";.*gene_name "(?P<gene_name>.*?)";').set_index('transcript_id')
        
        if quant in ["htseq", "featurecounts"]:
            exon = gtf.query('type == "exon"').copy()
            group = exon['attr'].str.extractall(r'gene_id "(?P<gene_id>.*?)";').reset_index(drop=True)
            group.index = exon.index
            exon['group'] = group
            exon.drop(columns=["type", "attr"], inplace=True)
            # all exon length
            gene_length = exon.groupby('group').apply(self.helper)
            gene_length = gene_length.to_frame('length')
            self.gene_info = pd.merge(self.id2gene, gene_length, left_index=True, right_index=True)
        
    def salmon_parser(self):
        data = []
        for file in self.input:
            df = pd.read_csv(file, sep="\t", index_col=0, usecols=['Name', 'EffectiveLength', 'TPM', 'NumReads'])
            name = file.split('/')[-2]
            df.columns = ["length", f"{name}_TPM", f"{name}_counts"]
            df[f'{name}_FPKM'] = df[f"{name}_counts"] * 10**9 / (df['length'] * df[f"{name}_counts"].sum())
            df.drop(columns=['length'], inplace=True)
            data.append(df)
        expression = pd.concat(data, axis=1)
        # transcript to gene level
        common = expression.index.intersection(self.trans2gene.index)
        self.exp = expression.groupby(self.trans2gene.loc[common, 'gene_name']).sum()

    def kallisto_parser(self):
        data = []
        for file in self.input:
            df = pd.read_csv(file, sep="\t", index_col=0, usecols=['target_id', 'eff_length', 'est_counts', 'tpm'])
            name = file.split('/')[-2]
            df.columns = ["length", f"{name}_counts", f"{name}_TPM"]
            df[f'{name}_FPKM'] = df[f"{name}_counts"] * 10**9 / (df['length'] * df[f"{name}_counts"].sum())
            df.drop(columns=['length'], inplace=True)
            data.append(df)
        expression = pd.concat(data, axis=1)
        # transcript to gene level
        common = expression.index.intersection(self.trans2gene.index)
        self.exp = expression.groupby(self.trans2gene.loc[common, 'gene_name']).sum()

    def rsem_parser(self):
        data = []
        for file in self.input:
            df = pd.read_csv(file, sep="\t", index_col=0, usecols=["gene_id", "expected_count", "TPM", "FPKM"])
            name = file.split('/')[-2]
            df.columns = [f"{name}_counts", f"{name}_TPM", f"{name}_FPKM"]
            data.append(df)
        expression = pd.concat(data, axis=1)
        # gene id to gene name
        common = expression.index.intersection(self.id2gene.index)
        self.exp = expression.loc[common]
        self.exp.index = self.id2gene.loc[common, 'gene_name']

    def featurecounts_parser(self):
        data = []
        for file in self.input:
            df = pd.read_csv(file, sep="\t", skiprows=1, index_col=0, usecols=[0, 6])
            df = pd.merge(self.gene_info, df, left_index=True, right_index=True).set_index('gene_name')
            name = file.split('/')[-2]
            df.columns = ["length", f"{name}_counts"]
            df[f'{name}_FPKM'] = df[f'{name}_counts'] * 10**9 / (df['length'] * df[f'{name}_counts'].sum())
            df[f'{name}_TPM'] = df[f'{name}_FPKM'] * 10**6 / df[f'{name}_FPKM'].sum()
            df.drop(columns=['length'], inplace=True)
            data.append(df)
        self.exp = pd.concat(data, axis=1)

    def htseq_parser(self):
        data = []
        for file in self.input:
            df = pd.read_csv(file, skipfooter=5, engine="python", sep="\t", header=None, names=["gene_id", "count"], index_col=0)
            df = pd.merge(self.gene_info, df, left_index=True, right_index=True).set_index('gene_name')
            name = file.split('/')[-2]
            df.columns = ["length", f"{name}_counts"]
            df[f'{name}_FPKM'] = df[f'{name}_counts'] * 10**9 / (df['length'] * df[f'{name}_counts'].sum())
            df[f'{name}_TPM'] = df[f'{name}_FPKM'] * 10**6 / df[f'{name}_FPKM'].sum()
            df.drop(columns=['length'], inplace=True)
            data.append(df)
        self.exp = pd.concat(data, axis=1)

    def run(self):
        self.logger.info(f"Write expression matrix to: {self.output}")
        self.exp.index.name = 'gene_name'
        index = self.exp.index.str.startswith('ENSG')
        lncRNA = self.exp.loc[index]
        prefix, suffix = os.path.splitext(self.output)
        if not lncRNA.empty:
            lncRNA.to_csv(prefix + ".lncRNA" + suffix, sep="\t")
        coding_gene = self.exp.loc[~index]
        coding_gene.to_csv(self.output, sep="\t")
        self.logger.info("Wlll done!")


if __name__ == "__main__":
    Wrapper(snakemake)