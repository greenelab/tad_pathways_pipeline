#!/bin/bash

# Example of a TAD_Pathways Analysis applied to Bone Mineral Density GWAS
set -o errexit

gwas_file='data/gwas_catalog/Bone_mineral_density_hg19.tsv'
snp_file='data/gwas_tad_snps/Bone_mineral_density_hg19_SNPs.tsv'
tad_file='data/gwas_tad_genes/Bone_mineral_density_hg19_SNPs_TAD_genelists.tsv'
trait='BMD'

# Perform WebGestalt pathway analysis and parse results
Rscript --vanilla scripts/webgestalt_run.R \
        --tad_genelist_file $tad_file \
        --output_name $trait

# Construct an evidence file - Nearest gene to gwas or not
python scripts/construct_evidence.py \
        --trait $trait \
        --gwas $gwas_file \
        --pathway 'skeletal system development'

# Summarize the evidence file
python scripts/summarize_evidence.py \
        --evidence 'results/BMD_gene_evidence.csv' \
        --snps $snp_file \
        --output_file 'results/BMD_gene_evidence_summary.tsv'

# Visualize overlap in TAD pathways curation
R --no-save --args 'results/BMD_gene_evidence.csv' \
        < scripts/integrative_summary.R 
