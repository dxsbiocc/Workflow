#!/bin/bash

# Setting root_dir in config.yaml
echo "Setting project dir ... "
script_path=$(dirname "$(realpath "$0")")
sed -i "s|{PLACEHOLDER}|$script_path|g" "$script_path/config/config.yaml"
echo "Done."

# Automatically find the Conda base directory
CONDA_BASE=$(conda info --base)

# Check if the path was successfully fetched
if [ -z "$CONDA_BASE" ]; then
    echo "Error: Unable to locate the Conda base directory."
    exit 1
fi

# Concatenate the full path of conda.sh
CONDA_SH="$CONDA_BASE/etc/profile.d/conda.sh"

# Initializing the Conda Environment
source $CONDA_SH

# Check if the snakemake virtual environment exists
env_exists=$(conda env list | grep "^snakemake\s")

if [ -z "$env_exists" ]; then
    # Install snakemake
    echo "Install snakemake environment ..."
    conda create -n snakemake python=3.10 snakemake=7.16 mamba -c bioconda -c conda-forge
else
    echo "Environment 'snakemake' exists."
fi

echo "Enter snakemake environment ..."
conda activate snakemake
# Install package
echo "Install packages: "
pip install -r ./requirement.txt
echo "Packages installed."

# cp scripts/base.py to snakemake-wrapper-utils
echo "Copy base.py ..."
location=$(pip show snakemake-wrapper-utils | awk '/Location/ {print substr($0, index($0,$2))}')
cp "$script_path/scripts/base.py" "$location/snakemake_wrapper_utils"
echo "Copy done."

echo "Exit snakemake environment"
conda deactivate