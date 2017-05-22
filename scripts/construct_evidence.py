"""
2016 Gregory Way
scripts/construct_evidence.py

Description:
Take as input the .tsv results from the WebGestalt analysis and outputs trait
specific TAD pathway genes

Usage:
Command line: python scripts/construct_evidence.py

    With the following required flags:

      --trait       a string of how the gestalt summary file was saved
      --gwas        the location of the GWAS nearest gene results
      --group            The group to subset the SNP file (optional)
      --pathway     a string of pathways to summarize (options are in the
                    gestalt summary files). To investigate multiple pathways
                    input a comma separated string.

                    For example:
                    --pathway 'skeletal system development'
                    --pathway 'skeletal system development,ossification'

Output:
a .csv evidence file (two columns: [gene, type of evidence])
"""

import os
import argparse
import pandas as pd

# Load Command Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-t', '--trait', help='symbol for trait data')
parser.add_argument("-g", "--gwas", help="location of gwas genelist")
parser.add_argument("-r", "--group", help="Group to subset evidence file",
                    default=None)
parser.add_argument("-p", "--pathway", help="pathway of interest")
args = parser.parse_args()

# Load Constants
trait = args.trait
gwas_file = args.gwas
gwas_group = args.group
try:
    pathway_df = pd.read_table(args.pathway)
    pathway = pathway_df[pathway_df['adjP'] < 0.05]['go_name'].tolist()
except ValueError:
    pathway = pathway.split(',')

trait_file = os.path.join('gestalt', '{}_complete_gestalt.tsv'.format(trait))
output_file = os.path.join('results', '{}_gene_evidence.csv'.format(trait))

# Load Data
pathway_genes = pd.read_csv(trait_file, delimiter='\t')
pathway_genes = pathway_genes.loc[pathway_genes['go_name'].isin(pathway), :]
pathway_genes = pathway_genes['symbol'].tolist()

gwas_genes_df = pd.read_table(gwas_file)
if gwas_group:
    gwas_genes_df = gwas_genes_df[gwas_genes_df['group'] == gwas_group]
gwas_genes = gwas_genes_df.MAPPED_GENE

# Process GWAS genes
full_gwas_genes = []
for gene in gwas_genes:
    gene = str(gene).replace(' - ', ', ').split(', ')
    full_gwas_genes += gene

# Process all genes
all_genes = set(full_gwas_genes + pathway_genes)
all_assignments = []

# Logic to determine genic evidence
for gene in all_genes:
    if gene in full_gwas_genes:
        assignment = 'gwas'
        if gene in pathway_genes:
            assignment = '{}_tad'.format(assignment)
    else:
        assignment = 'tad'
    all_assignments.append(assignment)

evidence = pd.DataFrame([all_genes, all_assignments]).T
evidence.columns = ['gene', 'evidence']
evidence.to_csv(output_file, sep=',', index=False)
