#!/bin/bash

# Example of a TAD_Pathways Analysis applied to Bone Mineral Density GWAS
set -o errexit

gwas_file='data/gwas_catalog/Bone_mineral_density_hg19.tsv'
snp_file='data/gwas_tad_snps/Bone_mineral_density_hg19_SNPs.tsv'
tad_file='data/gwas_tad_genes/Bone_mineral_density_hg19_SNPs_TAD_genelists.tsv'
trait='bmd'
evidence_file='results/bmd_gene_evidence.csv'

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
        --evidence $evidence_file \
        --snps $snp_file \
        --output_file 'results/BMD_gene_evidence_summary.tsv'

# Visualize overlap in TAD pathways curation
R --no-save --args $evidence_file \
        < scripts/integrative_summary.R 
