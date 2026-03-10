#!/bin/bash

# Get species and version
db=$1
rna="no"
spikein="no"

# Checking input parameters
if [ -z "$db" ]; then
    echo "Usage: $0 <species> [-RNA] [-spikein]"
    exit 1
fi

# Checking optional parameters
shift 2
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -RNA) rna="yes";;
        -spikein) spikein="yes";;
        hg19) ;;
        hg38) ;;
        mm10) ;;
        mm9) ;;
        *) echo "Unknown parameter: $1"; exit 1;;
    esac
    shift
done

# Create genome directory
if [ ! -d data/genome ]; then
    mkdir -p data/genome
fi

cd data/genome

fa=https://hgdownload.soe.ucsc.edu/goldenPath/$db/bigZips/$db.fa.gz
gtf=https://hgdownload.soe.ucsc.edu/goldenPath/$db/bigZips/genes/$db.refGene.gtf.gz

# Download genome data
echo "Downloading reference genome for $db ..."
if ! wget "$fa" -O "$db.fa.gz"; then
    echo "Error: Failed to download from $fa"
    exit 1
fi

echo "Downloading reference GTF for $db ..."
if ! wget "$gtf" -O "$db.refGene.gtf.gz"; then
    echo "Error: Failed to download from $gtf"
    exit 1
fi

if [[ "$rna" == "yes" ]]; then
    rna_fa=https://hgdownload.soe.ucsc.edu/goldenPath/$db/bigZips/refMrna.fa.gz

    echo "Downloading RNA reference genome for $db ..."
    if ! wget "$rna_fa" -O "$db.mrna.fa.gz"; then
        echo "Error: Failed to download from $rna_fa"
        exit 1
    fi
fi

# if [[ "$spikein" == "yes" ]]; then
#     spikein_fa=https://hgdownload.soe.ucsc.edu/goldenPath/$db/bigZips/chromFa.tar.gz
#     echo "Downloading spike-in (E.coli) reference genome ..."
#     if ! wget "$spikein_fa" -O "chromFa.tar.gz"; then
#         echo "Error: Failed to download from $spikein_fa"
#         exit 1
#     fi
# fi

ls ./*.gz | xargs gunzip