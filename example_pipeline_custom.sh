#!/bin/bash

# Example of a TAD_Pathways Analysis applied to a Custom SNP list
# For this example, the custom SNP list is the GWAS findings for
# Prostate Cancer. The data is used as a custom input.
set -o errexit

# Set constants
trait='custom'
custom_dir='results/custom_example/'
geneset='geneontology_Biological_Process'
pathway_name='skeletal_system_development'

# Set input file names
candidate_snp_file='custom_example.csv'

# Set output file names
candidate_snp_location_file=$custom_dir'custom_example_location.tsv'
candidate_snp_tad_file=$custom_dir'custom_example_tad_results.tsv'
nearest_gene_file=$custom_dir'custom_example_tad_results_nearest_gene.tsv'
pathway_file=$custom_dir'/'$trait'_gestalt.tsv'
evidence_file=$custom_dir'/'$trait'_gene_evidence.csv'
summary_file=$custom_dir'/'$trait'_gene_evidence_summary.tsv'

# Map SNPs to genomic location
Rscript --vanilla scripts/build_snp_list.R \
        --snp_file $candidate_snp_file \
        --output_file $candidate_snp_location_file

# Build a customized genelist based on SNP locations
python scripts/build_custom_tad_genelist.py \
        --snp_data_file $candidate_snp_location_file \
        --output_file $candidate_snp_tad_file

# Perform WebGestalt pathway analysis and parse results
Rscript --vanilla scripts/webgestalt_run.R \
        --tad_genelist_file $candidate_snp_tad_file \
        --output_name $trait \
        --output_directory $custom_dir

# Construct an evidence file - Nearest gene to gwas or not
python scripts/construct_evidence.py \
        --trait $trait \
        --gwas_file $nearest_gene_file \
        --pathway_file $pathway_file \
        --pathway "top" \
        --output_directory $custom_dir

# Summarize the evidence file
python scripts/summarize_evidence.py \
        --evidence $evidence_file \
        --snps $candidate_snp_tad_file \
        --output_file $summary_file

# Visualize overlap in TAD pathways curation
Rscript --vanilla scripts/integrative_summary.R \
        --evidence_file $evidence_file \
        --trait $trait \
        --output_directory $custom_dir
