#!/bin/bash

set -o errexit

# Example of a TAD_Pathways Analysis applied to Type 2 Diabetes GWAS

# After saving WebGestalt tsv file, parse its contents
python scripts/parse_gestalt.py --trait 't2d'

# Construct an evidence file - Nearest gene to gwas or not
python scripts/construct_evidence.py \
        --trait 't2d'\
        --gwas 'data/gwas_catalog/Type_2_diabetes_hg19.tsv'\
        --pathway 'peptide hormone secretion'

# Summarize the evidence file
python scripts/summarize_evidence.py \
        --evidence 'results/t2d_gene_evidence.csv' \
        --snps 'data/gwas_tad_snps/Type_2_diabetes_hg19_SNPs.tsv' \
        --output_file 'results/t2d_gene_evidence_summary.tsv'

# Visualize overlap in TAD pathways curation
R --no-save --args 'results/t2d_gene_evidence.csv' \
        < scripts/integrative_summary.R 
