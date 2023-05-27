import os
import re
import pandas as pd

def find_file(path):
    file_list = sorted(os.listdir(path))
    # fastq1
    p1 = re.compile('(.*)_R?1(_0[0-9]+)?.f[ast]{0,3}q.gz')
    # sample list
    sp_list = [p1.findall(f)[0][0] for f in file_list if p1.match(f)]
    #
    data = pd.DataFrame(sp_list, columns=['sample'])
    data['fastq1'] = [f for f in file_list if p1.match(f)]
    # fastq2
    p2 = re.compile('(.*)_R?2(_0[0-9]+)?.f[ast]{0,3}q.gz')
    fq2_list = [f for f in file_list if p2.match(f)]
    if len(fq2_list) == 0:
        data['fastq2'] = None
    elif len(sp_list) != len(fq2_list):
        raise ValueError('The number of files in fq1 and fq2 is not the same!')
    else:
        data['fastq2'] = fq2_list
    return data