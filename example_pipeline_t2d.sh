#!/bin/bash

# Example of a TAD_Pathways Analysis applied to Type 2 Diabetes GWAS
set -o errexit

gwas_file='data/gwas_catalog/Type_2_diabetes_hg19.tsv'
snp_file='data/gwas_tad_snps/Type_2_diabetes_hg19_SNPs.tsv'
tad_file='data/gwas_tad_genes/Type_2_diabetes_hg19_SNPs_TAD_genelists.tsv'
trait='t2d'
evidence_file='results/t2d_gene_evidence.csv'
pathway_p_values_file='gestalt/t2d_pvals.tsv'

# Perform WebGestalt pathway analysis and parse results
Rscript --vanilla scripts/webgestalt_run.R \
        --tad_genelist_file $tad_file \
        --output_name $trait

# Construct an evidence file
# In this case, select all significant pathways, no preselection required
python scripts/construct_evidence.py \
        --trait $trait \
        --gwas $gwas_file \
        --pathway $pathway_p_values_file

# Summarize the evidence file
python scripts/summarize_evidence.py \
        --evidence $evidence_file \
        --snps $snp_file \
        --output_file 'results/t2d_gene_evidence_summary.tsv'

# Visualize overlap in TAD pathways curation
R --no-save --args $evidence_file \
        < scripts/integrative_summary.R 
