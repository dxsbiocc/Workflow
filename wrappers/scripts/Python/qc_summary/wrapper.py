# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2024/01/20 15:41:36
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : Summary QC information
# ============================================================


import re
import os
from glob import glob

import pandas as pd
from snakemake import shell
from snakemake_wrapper_utils.base import WrapperBase, get_logger


class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        self.logger = get_logger(snakemake.rule, filename=str(snakemake.log))
        super().__init__(snakemake)

    def parse_mapping(self, file_list):
        mapping_data, spikein = [], []
        for file in file_list:
            name = re.findall("bowtie2_(.+?)(_spikein)?.summary", file)[0][0]
            with open(file) as fp:
                for line in fp:
                    if 'reads; of these:' in line:
                        total = int(re.findall(r'(\d+) reads; of these:', line)[0])
                    if 'exactly 1 time' in line:
                        uniq_align = int(re.findall(r'\s+(\d+) \(.+\) aligned concordantly exactly 1 time', line)[0])
                    if ' >1 times' in line:
                        multi_align = int(re.findall(r'\s+(\d+) \(.+\) aligned concordantly >1 times', line)[0])
                mapping_ratio = '{:.2%}'.format((uniq_align+multi_align)/total)
                if 'spikein.summary' in file:
                    spikein.append([name, mapping_ratio])
                else:
                    mapping_data.append([name, total * self.factor, uniq_align * self.factor, mapping_ratio])
        mapping_data = pd.DataFrame(mapping_data, columns=['sample', 'total_reads', 'unique_align', 'mapping_ratio']).set_index('sample')
        spikein = pd.DataFrame(spikein, columns=['sample', 'spikein_ratio']).set_index('sample')
        mapping_data['spikein_ratio'] = spikein
        return mapping_data

    def parse_dedup_metric(self, file_list):
        metric = []
        for file in file_list:
            name = file.split('/')[-2]
            df = pd.read_csv(file, sep='\t', skiprows=6, nrows=1).squeeze()
            dup_ratio = '{:.2%}'.format(df['PERCENT_DUPLICATION'])
            size = df['ESTIMATED_LIBRARY_SIZE']
            metric.append([name, dup_ratio, size])
        return pd.DataFrame(metric, columns=['sample', 'duplication_ratio', 'estimated_library_size']).set_index('sample')
    
    def parse_FRIP(self, file_list):
        frip = []
        for file in file_list:
            name = os.path.basename(file).replace('.txt', '')
            df = pd.read_csv(file, sep='\t')
            df.index = [name]
            frip.append(df)
        frip = pd.concat(frip)
        frip.columns = ['used_reads', 'peak_reads', 'FRiP']
        return frip
    
    def parse_region_frac(self, file_list):
        region = []
        for file in file_list:
            name = os.path.basename(file).replace('.txt', '')
            df = pd.read_csv(file, sep='\t')
            df.index = [name]
            region.append(df)
        region = pd.concat(region)
        region.columns = ['reads_all', 'reads_in_dnase', 'reads_in_blacklist', 'reads_in_promoter', 'reads_in_enhancer']
        return region

    def parse_library_complexity(self, file_list):
        libs = []
        for file in file_list:
            name = os.path.basename(file).replace('.txt', '')
            df = pd.read_csv(file, sep='\t')
            df.index = [name]
            libs.append(df)
        libs = pd.concat(libs)
        libs.columns = ['total_fragments', 'distinct_fragments', 'read1', 'read2', 'NRF', 'PBC1', 'PBC2']
        return libs

    def parse_mito_frac(self, file_list):
        idxstats = []
        for file in file_list:
            name = os.path.basename(file).replace('.idxstats', '')
            df = pd.read_csv(file, sep='\t', usecols=[0, 2], header=None, names=['chrom', 'mitochondria'])
            tmp = df.query('chrom == "chrM"').drop(columns='chrom') / df['mitochondria'].sum()
            tmp.index = [name]
            idxstats.append(tmp)
        idxstats = pd.concat(idxstats)
        return idxstats

    def parse_peak_nano(self, file_list):
        anno = []
        for file in file_list:
            name = os.path.basename(file).replace('.peakAnno.txt', '')
            df = pd.read_csv(file, sep='\t', usecols=['annotation'])
            a = df.annotation.str.extractall(r"(?P<region>\w+) \(.*\)").value_counts()
            b = df.loc[~df.annotation.str.contains(r'\(')].value_counts()
            df = pd.concat([a, b])
            df.name = name
            anno.append(df)
        anno = pd.concat(anno, axis=1).transpose()
        anno.columns = [f"peak_in_{a[0].lower()}" for a in anno.columns]
        return anno

    def parser(self):
        self.output = str(self.snakemake.output)
        # paired factor
        paired = self.snakemake.params.get('paired')
        self.factor = 2 if paired else 1
        # output data
        data = []
        # mapping
        mapping = self.snakemake.params.get('mapping')
        if mapping == 'bowtie2':
            self.logger.info("parse bowtie2 mapping result ...")
            mapping_files = self.snakemake.input.get('mapping_files')
            data.append(self.parse_mapping(mapping_files))
            self.logger.info("parse successed")
        else:
            self.logger.error(f'unsupported mapping tool: {mapping}')
        # dedup
        dedup = self.snakemake.params.get('dedup')
        if dedup == 'markduplicates':
            self.logger.info("parse markduplicated result ...")
            metric_files = self.snakemake.input.get('dedup_fils')
            data.append(self.parse_dedup_metric(metric_files))
            self.logger.info("parse successed")
        else: 
            self.logger.error(f'unsupported dedup tool: {dedup}')
        # FRiP
        self.logger.info("parse FRiP files ...")
        frip_files = self.snakemake.input.get('frip_files')
        data.append(self.parse_FRIP(frip_files))
        self.logger.info("parse successed")
        # annotate regions fraction
        self.logger.info("parse annotated regions ...")
        regiond_files = self.snakemake.input.get('isize_files')
        data.append(self.parse_region_frac(regiond_files))
        self.logger.info("parse successed")
        # library complexity
        self.logger.info("parse library complexity files ...")
        lib_complexity_files = self.snakemake.input.get('lib_complexity_files')
        data.append(self.parse_library_complexity(lib_complexity_files))
        self.logger.info("parse successed")
        # mitochondria percent
        self.logger.info("parse mapping idxstats files ...")
        idxstats_files = self.snakemake.input.get('idxstats_files')
        data.append(self.parse_mito_frac(idxstats_files))
        self.logger.info("parse successed")
        # anno
        self.logger.info("parse peak annotation files ...")
        anno_files = self.snakemake.input.get('anno_files')
        data.append(self.parse_peak_nano(anno_files))
        self.logger.info("parse successed")

        self.data = pd.concat(data, axis=1)
        self.data.index.name = 'sample'
    
    def run(self):
        if not os.path.exists(self.output):
            self.logger.info('touch file')
            shell("touch {self.output}")
        self.logger.info(f"writing result to {self.output}")
        if self.output.endswith('xlsx'):
            self.data.to_excel(self.output)
        elif self.output.endswith('csv'):
            self.data.to_csv(self.output)
        elif self.output.endswith('txt') or self.output.endswith('xls') or self.output.endswith('tsv'):
            self.data.to_csv(self.output, sep='\t')
        else:
            self.logger.error('unsupport file format, using default format(.xlsx)')
            self.data.to_excel(self.output + ".xlsx")
        self.logger.info('Well done!')
        

if __name__ == '__main__':
    Wrapper(snakemake)