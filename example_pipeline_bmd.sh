#!/bin/bash

set -o errexit

# Example of a TAD_Pathways Analysis applied to Bone Mineral Density GWAS

# After saving WebGestalt tsv file, parse its contents
python scripts/parse_gestalt.py --trait 'bmd'

# Construct an evidence file - Nearest gene to gwas or not
python scripts/construct_evidence.py \
        --trait 'bmd'\
        --gwas 'data/gwas_catalog/Bone_mineral_density_hg19.tsv'\
        --pathway 'skeletal system development'

# Summarize the evidence file
python scripts/summarize_evidence.py \
        --evidence 'results/bmd_gene_evidence.csv' \
        --snps 'data/gwas_tad_snps/Bone_mineral_density_hg19_SNPs.tsv' \
        --output_file 'results/bmd_gene_evidence_summary.tsv'

# Visualize overlap in TAD pathways curation
R --no-save --args 'results/bmd_gene_evidence.csv' \
        < scripts/integrative_summary.R 
