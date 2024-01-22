# -*- encoding: utf-8 -*-
# ============================================================
# File        : wrapper.py
# Time        : 2024/01/18 09:33:03
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : The role of the current file 
# ============================================================


import os
import warnings
import numpy as np
import pandas as pd
from scipy.signal import find_peaks_cwt
from snakemake_wrapper_utils.base import WrapperBase, get_logger
from plotnine import (ggplot, aes, geom_col, scale_x_continuous, 
                      labs, theme_classic, theme, element_text)

warnings.filterwarnings("ignore")

class Wrapper(WrapperBase):

    def __init__(self, snakemake) -> None:
        self.logger = get_logger(snakemake.rule, filename=str(snakemake.log))
        super().__init__(snakemake)

    def parser(self):
        infile = str(self.snakemake.input)
        self.name = os.path.basename(infile).replace('.hist_data.txt', '')
        self.data = pd.read_csv(infile, sep='\t', skiprows=11, names=["isize", "fr_count"])
    
    def fragment_length_qc(self):
        results = []

        NFR_UPPER_LIMIT = 150
        MONO_NUC_LOWER_LIMIT = 150
        MONO_NUC_UPPER_LIMIT = 300

        data = self.data.values
        # % of NFR vs res
        nfr_reads = data[data[:, 0] < NFR_UPPER_LIMIT][:, 1]
        percent_nfr = round(nfr_reads.sum() / data[:, 1].sum(), 4)
        results.append(["Fraction of reads in NFR", percent_nfr >= 0.4, 
                        "OK" if percent_nfr >= 0.4 else f"{percent_nfr} out of range [0.4, INF]"])

        # % of NFR vs mononucleosome
        mono_nuc_reads = data[
            (data[:, 0] > MONO_NUC_LOWER_LIMIT) &
            (data[:, 0] <= MONO_NUC_UPPER_LIMIT)][:, 1]

        percent_nfr_vs_mono_nuc = round(nfr_reads.sum() / mono_nuc_reads.sum(), 4)
        results.append(["NFR / mono-nuc reads", percent_nfr_vs_mono_nuc >= 2.5, 
                        "OK" if percent_nfr_vs_mono_nuc >= 2.5 else f"{percent_nfr_vs_mono_nuc} out of range [2.5, INF]"])

        # peak locations
        pos_start_val = data[0, 0]  # this may be greater than 0
        peaks = find_peaks_cwt(data[:, 1], np.array([25]))
        nuc_range_metrics = [
            ('Presence of NFR peak', 20 - pos_start_val, 90 - pos_start_val),
            ('Presence of Mono-Nuc peak',
            120 - pos_start_val, 250 - pos_start_val),
            ('Presence of Di-Nuc peak',
            300 - pos_start_val, 500 - pos_start_val)]
        for range_metric in nuc_range_metrics:
            if len([peak for peak in peaks if range_metric[1] <= peak <= range_metric[2]]):
                results.append([range_metric[0], True, "OK"])
            else:
                results.append([range_metric[0], False, f"Cannot find element in range [{range_metric[1]}, {range_metric[2]}]"])
        
        return results
    
    def plot(self):
        plot = (ggplot(self.data) 
             + geom_col(aes("isize", "fr_count"), colour = "#60c5ba") 
             + scale_x_continuous(expand = (0, 0, 0, 0), limits = (0, 1000)) 
             + labs(x = "Insert Size", y = "Counts", title = self.name) 
             + theme_classic() 
             + theme(plot_title = element_text(hjust = 0.5)))
        return plot
    
    def run(self):
        # qc file
        self.logger.info("fragment length qc")
        results = self.fragment_length_qc()
        outfile = self.snakemake.output.get('qc')
        df = pd.DataFrame(results, columns=['metric', 'qc_pass', 'message'])
        df.to_csv(outfile, sep='\t', index=False)
        # plot
        self.logger.info("plot fragment length")
        pdf = self.snakemake.output.get('pdf')
        p = self.plot()
        p.save(pdf, format="pdf", height=6, width=8)
        self.logger.info("Well done")
    

if __name__ == '__main__':
    Wrapper(snakemake)