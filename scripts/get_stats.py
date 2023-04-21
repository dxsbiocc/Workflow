import re
import os
import pandas as pd


def get_mapping(file, spikein=False):
    with open(file) as fp:
        for line in fp:
            if 'reads; of these:' in line:
                total = int(re.findall('(\d+) reads; of these:', line)[0])
            if 'exactly 1 time' in line:
                uniq_align = int(re.findall('\s+(\d+) \(.+\) aligned concordantly exactly 1 time', line)[0])
            if ' >1 times' in line:
                multi_align = int(re.findall('\s+(\d+) \(.+\) aligned concordantly >1 times', line)[0])
        mapping_ratio = '{:.2%}'.format((uniq_align+multi_align)/total)
        if spikein:
            mapping = [uniq_align+multi_align, mapping_ratio]
        else:
            mapping = [total, uniq_align, mapping_ratio]
    return mapping

def get_metric(file):
    metric = pd.read_csv(file, sep='\t', skiprows=6, nrows=1).squeeze()
    dup_ratio = '{:.2%}'.format(metric['PERCENT_DUPLICATION'])
    size = metric['ESTIMATED_LIBRARY_SIZE']
    return [dup_ratio, size]

def get_FRiP(bam, peak):
    mapped_reads =  int(os.popen(f'samtools view -c {bam}').read().strip())
    reads_on_peak = int(os.popen(f"bedtools sort -i {peak} \
                              | bedtools merge -i stdin | bedtools intersect -u -nonamecheck \
                              -a {bam} -b stdin -ubam | samtools view -c").read().strip())
    frip = '{:.2%}'.format(reads_on_peak/mapped_reads)
    return [mapped_reads, frip]


def main(snakemake):
    columns = [
        'Sample', 'TotalReads', 'UniqueMapping', 'MappingRate', 'MappedSpikein', 
        'MappingRateSpikein', 'DuplicationRate', 'LibrarySize', 'FilteredReads', 'FRiP'
    ]

    frip = {}
    for files in zip(*[sorted(files) for files in [snakemake.input.bam, snakemake.input.peak]]):
        pair = files[0].split('/')[-1].split('.')[0]
        frip[pair] = get_FRiP(files[0], files[1])

    res = []
    file_list = [sorted(files) for files in [snakemake.params.sample_list, snakemake.input.summary, 
                                                snakemake.input.spikein, snakemake.input.metric]]
    for files in zip(*file_list):
        record = [files[0]]
        record.extend(get_mapping(files[1]))
        record.extend(get_mapping(files[2], spikein=True))
        record.extend(get_metric(files[3]))
        if (files[0] in frip):
            record.extend(frip[files[0]])
        else:
            record.extend(['-', '-'])
        res.append(record)

    data = pd.DataFrame(res, columns=columns)
    data.to_csv(snakemake.output[0], index=False)


if __name__ == '__main__':
    main(snakemake)