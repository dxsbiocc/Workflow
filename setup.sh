#!/bin/bash

# setting root_dir in config.yaml
echo "Setting project dir ... "
script_path=$(dirname "$(realpath "$0")")
sed -i "s|{PLACEHOLDER}|$script_path|g" "$script_path/config/config.yaml"
echo "Done."

# install snakemake
echo "Install snakemake environment ..."
conda create -n snakemake python=3.10 snakemake=7.16 mamba -c bioconda -c conda-forge

echo "Enter snakemake environment ..."
conda activate snakemake
# install package
echo "Install packages: "
pip install -r ./requirement.txt
echo "Packages installed."

# cp utils/base.py to snakemake-wrapper-utils
echo "Copy base.py ..."
location=$(pip show snakemake-wrapper-utils | awk '/Location/ {print substr($0, index($0,$2))}')
cp "$script_path/utils/base.py" "$location/snakemake_wrapper_utils"
echo "Copy done."

echo "Exit snakemake environment"
conda deactivate