#!/bin/bash

# install package
echo "\e[32m install packages: \e[0m"
pip install -r ./requirement.txt
echo "\e[32m Done \e[0m"

# cp utils/base.py to snakemake-wrapper-utils
echo "\e[32m Copy base.py \e[0m"
location=$(pip show snakemake-wrapper-utils | awk '/Location/ {print substr($0, index($0,$2))}')
cp utils/base.py "$location"
echo "\e[32m Done \e[0m"

# Check if the bedtools is installed
if ! command -v mamba &> /dev/null; then
    echo "\e[32m if 'mamba' does not exist, will use conda install... \e[0m"
    conda install -n base -c conda-forge mamba
else
    echo "\e[32m mamba already installed! \e[0m"
fi

# Check if the bedtools is installed
if ! command -v bedtools &> /dev/null; then
    echo "\e[32m if 'bedtools' does not exist, will use conda install... \e[0m"
    conda install -c bioconda bedtools
else
    echo "\e[32m bedtools already installed! \e[0m"
fi

# setting root_dir in config.yaml
echo "\e[32m Setting project dir ... \e[0m"
script_path=$(dirname "$(realpath "$0")")
sed -i "s|{PLACEHOLDER}|$script_path|g" "$script_path/config/config.yaml"
echo "\e[32m Done \e[0m"
