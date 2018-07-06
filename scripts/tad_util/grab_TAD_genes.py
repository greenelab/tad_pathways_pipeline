"""
2016 Gregory Way
TAD Pathways
scripts/tad_util/grab_TAD_genes.py

Description:
Uses intemediate files and gene assignment in TADs index file to curate
a final list of all trait-specific genes that fall in signal TADs

Usage:
Is called by 'scripts/run_pipeline.sh'

Output:
Trait specific .tsv files of one column each indicating all the genes that fall
in signal TADs
"""

import argparse
import pandas as pd

# Load Command Arguments

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--tad_gwas_file", help="Location of TAD/GWAS data")
parser.add_argument("-o", "--output_tad_file",
                    help="Name of the full genelist output file")
parser.add_argument("-g", "--output_gwas_file",
                    help="Name of the GWAS genelist output file")
args = parser.parse_args()

# Load Constants
tad_gwas_loc = args.tad_gwas_file
output_file = args.output_tad_file

# Load data
tad_genes_df = pd.read_table('index/GENE_index_hg19_hESC.tsv.bz2', index_col=0)
tad_gwas_df = pd.read_table(tad_gwas_loc)
tad_gwas_df = tad_gwas_df[~tad_gwas_df.TADidx.isnull()]
tad_gwas_df.TADidx = tad_gwas_df.TADidx.astype(int)

# Get all TAD based genes
tad_ids = map(str, tad_gwas_df['TADidx'].tolist())
tad_genes_df = tad_genes_df[tad_genes_df["TAD_id"].isin(tad_ids)]\
                                                  .sort_values(by='gene_name')

tad_genes_df.to_csv(output_file, index=False, header=True, sep='\t')
