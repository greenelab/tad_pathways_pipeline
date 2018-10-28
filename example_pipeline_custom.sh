#!/bin/bash

# Example of a TAD_Pathways Analysis applied to a Custom SNP list
# For this example, the custom SNP list is the GWAS findings for
# Prostate Cancer. The data is used as a custom input.
set -o errexit

# Define filenames
candidate_snp_file='custom_example.csv'
candidate_snp_location_file='results/custom_example_location.tsv'
candidate_snp_tad_file='results/custom_example_tad_results.tsv'
nearest_gene_file='results/custom_example_tad_results_nearest_gene.tsv'
trait='custom'
evidence_file='results/custom_gene_evidence.csv'
pathway_file='pathway_results/custom_enrichr_results.tsv'

# Map SNPs to genomic location
Rscript --vanilla scripts/build_snp_list.R \
        --snp_file $candidate_snp_file \
        --output_file $candidate_snp_location_file

# Build a customized genelist based on SNP locations
python scripts/build_custom_tad_genelist.py \
        --snp_data_file $candidate_snp_location_file \
        --output_file $candidate_snp_tad_file

# Perform pathway analysis using gseapy
python scripts/pathway_analysis.py \
        --tad_genelist_file $candidate_snp_tad_file \
        --output_name $trait \
        --output_directory "pathway_results" \
        --gene_sets "KEGG_2016"

# Construct an evidence file - Nearest gene to gwas or not
python scripts/construct_evidence.py \
            --trait $trait \
            --gwas_file $nearest_gene_file \
            --pathway_file $pathway_file \
            --all_sig_pathways

# Summarize the evidence file
python scripts/summarize_evidence.py \
            --evidence $evidence_file \
            --snps $candidate_snp_tad_file \
            --output_file 'results/custom_gene_evidence_summary.tsv'

# Visualize overlap in TAD pathways curation
R --no-save --args $evidence_file \
                < scripts/integrative_summary.R
