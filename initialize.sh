#!/bin/bash

set -o errexit

# Download necessary index and curation files

url=https://zenodo.org/record/1311211/files/tad_pathways_data.tar.gz
wget $url

# Extract data
tar -zxvf tad_pathways_data.tar.gz

# Stabilize computing environments in Python and R
conda env create --quiet --force --file environment.yml
conda activate tad_pathways

Rscript scripts/install.R
