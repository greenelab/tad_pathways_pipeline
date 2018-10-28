"""
2018 Gregory Way
scripts/pathway_analysis.py

Description:
Runs pathway analysis on Input genelist

Usage:
The script is run in the command line

           python scripts/pathway_analysis.py

    with command line arguments:

           --tad_genelist_file       <LOCATION_OF_GWAS_FILE.csv>
           --output_name             <FILENAME_ABBREVIATION>
           --output_directory        <DIRECTORY TO SAVE RESULTS>
           --gene_sets               <STRING IDENTIFIERS OF GENE LISTS>

     An example gene_set input is "KEGG_2016". For a full list view:
     http://amp.pharm.mssm.edu/Enrichr/#stats

Output:
Overrepresented pathways represented in the given input SNP list
"""

import os
import argparse
import pandas as pd
import gseapy as gp

# Load Command Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--tad_genelist_file",
                    help="Location of TAD genes")
parser.add_argument("-o", "--output_name",
                    help="Abbrev. to save results"),
parser.add_argument("-d", "--output_directory",
                      help="Directory to save results",
                      default="pathway_results")
parser.add_argument("-g", "--gene_sets", nargs='+',
                      help="Geneset names to perform enrichment analysis over",
                      default="KEGG_2016")
args = parser.parse_args()

# Load arguments
tad_gene_file = args.tad_genelist_file
output_name = args.output_name
output_dir = args.output_directory
gene_sets = args.gene_sets

# Save output file locations
out_dir = os.path.join(output_dir, output_name)

# Load data
gene_df = pd.read_table(tad_gene_file)
genes = gene_df.gene_name
output_results_file = "{}_enrichr_results.tsv".format(out_dir)

# Perform the enrichment analysis
enr = gp.enrichr(gene_list=genes,
                 description=output_name,
                 gene_sets=gene_sets,
                 outdir=out_dir,
                 cutoff=0.5)

# Output results
results_df = enr.results.sort_values(by='Combined Score', ascending=False)
results_df.to_csv(output_results_file, sep='\t', index=False)
