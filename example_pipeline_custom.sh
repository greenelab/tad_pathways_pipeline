#!/bin/bash

# Example of a TAD_Pathways Analysis applied to a Custom SNP list
# For this example, the custom SNP list is the GWAS findings for
# Prostate Cancer. The the data is used as a custom input.

# Map SNPs to genomic location
Rscript --vanilla scripts/build_snp_list.R \
        --snp_file 'custom_example.csv' \
        --output_file 'results/custom_example_location.tsv'

# Build a customized genelist to input into WebGestalt
python scripts/build_custom_tad_genelist.py \
        --snp_data_file 'results/custom_example_location.tsv' \
        --output_file 'results/custom_example_tad_results.tsv'

# After saving WebGestalt tsv file, parse its contents
python scripts/parse_gestalt.py --trait 'custom'

# Construct an evidence file - Nearest gene to gwas or not
python scripts/construct_evidence.py \
            --trait 'custom'\
            --gwas 'results/custom_example_tad_results_nearest_gene.tsv'\
            --pathway 'epidermis development,antigen processing and presentation'

# Summarize the evidence file
python scripts/summarize_evidence.py \
            --evidence 'results/custom_gene_evidence.csv' \
            --snps 'results/custom_example_tad_results.tsv' \
            --output_file 'results/custom_gene_evidence_summary.tsv'

# Visualize overlap in TAD pathways curation
R --no-save --args 'results/custom_gene_evidence.csv' \
                < scripts/integrative_summary.R 

