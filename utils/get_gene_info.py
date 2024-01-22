# -*- encoding: utf-8 -*-
# ============================================================
# File        : get_gene_info.py
# Time        : 2024/01/19 15:47:18
# Author      : dengxsh
# Version     : 1.0
# Contact     : 920466915@qq.com
# License     : MIT
# Copyright   : Copyright (c) 2022, dengxsh
# Description : 
#   This script is designed to extract gene information from a 
#   GTF (Gene Transfer Format) file. 
# ============================================================


import os
import argparse
import pandas as pd

from utils import get_logger, close_logger

logger = get_logger(name="get_gene_info")

def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC filter.')
    parser.add_argument("-g", "--gtf", dest="gtf", type=str,
                        help='GTF file')
    parser.add_argument("-t", "--type", dest="type", type=str, choices=['all', 'gene', 'transcript'],
                        default='gene', help="gene type.")
    parser.add_argument("-a", "--attrs", dest='attrs', nargs='*', default=["gene_name", "hgnc_id"],
                        help="attributes of gene")

    args = parser.parse_args()

    logger.info(args)
    return args

def get_info(gtf, type='all', attrs=[]):
    """
    Get gene attributes information from GTF file, such as gene location(chrom, start, end)

    Parameters
    ==========
    gtf : DataFrame object, read from GTF file
    type : 
        * all: 
        * gene:
        * transcript
    attrs : str or list, or array-like
        * gene_id:
        * gene_name:
    """
    pass

def main():
    args = parse_arguments()
    if not os.path.exists(args.gtf):
        # raise FileNotFoundError(f'{args.gtf} not found!')
        logger.error(f'{args.gtf} not found!')
        # coler handler
        close_logger(logger)
        exit(1)
    logger.info(f'Read GTF file: {args.gtf}')
    gtf = pd.read_csv(args.gtf, skiprows=5, sep="\t", usecols=[0, 2, 3, 4, 8], header=None, 
                      names=['chrom', 'type', 'start', 'end', 'attr'])
    
    # coler handler
    close_logger(logger)
    


if __name__ == '__main__':
    main()