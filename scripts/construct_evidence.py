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
      --group       The group to subset the SNP file (optional)
      --pathway     a string of pathways to summarize (options are in the
                    gestalt summary files). To investigate multiple pathways
                    input a comma separated string.

                    For example:
                    --pathway 'skeletal system development'
                    --pathway 'skeletal system development,ossification'

    And the following optional flags:

      --results_directory       the directory of where to save results
      --gestalt_directory       the directory to load gestalt results from
      --all_pathways            will subset to all significant pathways
                                and not just consider the top pathway

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
parser.add_argument("-d", "--results_directory", help="where to save results",
                    default="results")
parser.add_argument("-e", "--gestalt_directory",
                    help="location of gestalt results", default="gestalt")
parser.add_argument("-a", "--all_pathways", action='store_true',
                    help="consider all pathways in evidence")
args = parser.parse_args()

# Load Constants
trait = args.trait
gwas_file = args.gwas
gwas_group = args.group
all_pathways = args.all_pathways
try:
    pathway_df = pd.read_table(args.pathway)
    pathway = pathway_df.query('adjP < 0.05').loc[:, 'go_name'].tolist()
except:
    pathway = args.pathway.split(',')

results_dir = args.results_directory
gestalt_dir = args.gestalt_directory

trait_file = os.path.join(gestalt_dir, '{}_gestalt.tsv'.format(trait))
if all_pathways:
    trait = '{}_all-sig-pathways'.format(trait)
output_file = os.path.join(results_dir, '{}_gene_evidence.csv'.format(trait))

# Load Data
drop_col = ['link', 'count', 'observed', 'expected', 'R', 'overlapGene']
pathway_genes_df = (
    pd.read_csv(trait_file, delimiter='\t')
    .query('go_name in @pathway')
    .drop(drop_col, axis='columns')
    )

pathway_genes = pathway_genes_df.loc[:, 'symbol'].tolist()

gwas_genes_df = pd.read_table(gwas_file)
if gwas_group:
    gwas_genes_df = gwas_genes_df.query('group == @gwas_group')
gwas_genes = gwas_genes_df.MAPPED_GENE

# Process GWAS genes
full_gwas_genes = []
for gene in gwas_genes:
    gene = str(gene).replace(' - ', ', ').split(', ')
    full_gwas_genes += gene

# Process all genes
all_genes = set(full_gwas_genes + pathway_genes)

# Logic to determine genic evidence
all_assignments = []
for gene in all_genes:
    if gene in full_gwas_genes:
        assignment = 'gwas'
        if gene in pathway_genes:
            assignment = '{}_tad'.format(assignment)
    else:
        assignment = 'tad'
    all_assignments.append(assignment)

evidence_df = pd.DataFrame([all_genes, all_assignments]).T
evidence_df.columns = ['gene', 'evidence']

# Write output to file
evidence_df = evidence_df.merge(pathway_genes_df, how='left', left_on='gene',
                                right_on='symbol')

if not all_pathways:
    evidence_df = (
        evidence_df
        .sort_values(by='pval')
        .groupby('gene', as_index=False)
        .first()
        )

(
    evidence_df
    .merge(gwas_genes_df, how='left', left_on='gene', right_on='MAPPED_GENE')
    .set_index('gene')
    .drop(['Unnamed: 0', 'MAPPED_GENE', 'symbol'], axis='columns')
    .to_csv(output_file, sep=',', index=True)
)
