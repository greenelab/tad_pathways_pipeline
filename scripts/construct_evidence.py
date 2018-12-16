"""
Gregory Way 2018
scripts/construct_evidence.py

Description:
Take as input the .tsv results from the WebGestalt analysis and outputs trait
specific TAD pathway genes

Usage:
Command line: python scripts/construct_evidence.py

    With the following required flags:

      --trait          a string of how the gestalt summary file was saved
      --gwas_file      the location of the GWAS nearest gene results
      --pathway_file   the location of the results of the GSEA analysis

    And the following optional flags:

      --pathway                 a string of pathways to summarize.
                                To investigate multiple pathways input a comma
                                separated string. Defaults to "top".
                        For example:
                        --pathway 'skeletal system development'
                        --pathway 'skeletal system development,ossification'
      --gwas_group              The group to subset the SNP file
      --all_sig_pathways        will subset to all significant pathways
                                and not just consider the top pathway
      --pathway_sig_cutoff      the alpha value to consider pathway significant
      --output_directory        the directory of where to save results

Output:
a .csv evidence file determining if the gene is located in the TAD or is the
nearest gene in a GWAS analysis
"""

import os
import argparse
import pandas as pd

# Load Command Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--trait", help="symbol for trait data")
parser.add_argument("-g", "--gwas_file", help="location of gwas genelist")
parser.add_argument("-r", "--gwas_group", help="Group to subset evidence file",
                    default=None)
parser.add_argument("-f", "--pathway_file",
                    help="file storing results of the pathway analysis")
parser.add_argument("-p", "--pathway", help="pathway of interest", default="top")
parser.add_argument("-a", "--all_sig_pathways", action="store_true",
                    help="consider all significant pathways in evidence")
parser.add_argument("-c", "--pathway_sig_cutoff", default=0.05,
                    help="Adjusted p-value cutoff for enrichment significance")
parser.add_argument("-o", "--output_directory", help="where to save results",
                    default="results")
args = parser.parse_args()

# Load Constants
trait = args.trait
gwas_file = args.gwas_file
gwas_group = args.gwas_group
pathway_file = args.pathway_file
pathway = args.pathway
all_sig_pathways = args.pathway
cutoff = args.pathway_sig_cutoff
output_dir = args.output_directory

# Load pathway analysis results data and subset based on target pathway(s)
pathway_results_df = pd.read_table(pathway_file)
if all_sig_pathways:
    pathway_results_df = pathway_results_df[
        pathway_results_df["adjP"] < cutoff
    ]
    trait = "{}_all-sig-pathways".format(trait)
elif pathway == "top":
    pathway_results_df = pathway_results_df.sort_values(by='adjP', ascending=True)
    top_pathway = pathway_results_df.id.tolist()[0]
    pathway_results_df = pathway_results_df.query("id == @top_pathway")
    trait = "{}_top-pathway".format(trait)
elif pathway == "all":
    next
else:
    pathway_results_df = pathway_results_df.query("term == @pathway")
    trait = "{}_{}".format(trait, pathway)

# Setup ouput files
output_file = os.path.join(output_dir, "{}_gene_evidence.csv".format(trait))

# Get all the pathway genes of interest
pathway_genes = pathway_results_df.symbol.tolist()

gwas_genes_df = pd.read_table(gwas_file)
if gwas_group:
    gwas_genes_df = gwas_genes_df.query("group == @gwas_group")
gwas_genes = gwas_genes_df.MAPPED_GENE

# Process GWAS genes
full_gwas_genes = []
for gene in gwas_genes:
    gene = str(gene).replace(" - ", ", ").split(", ")
    full_gwas_genes += gene

# Process all genes
all_genes = set(full_gwas_genes + pathway_genes)

# Logic to determine genic evidence
all_assignments = []
for gene in all_genes:
    if gene in full_gwas_genes:
        assignment = "gwas"
        if gene in pathway_genes:
            assignment = "{}_tad".format(assignment)
    else:
        assignment = "tad"
    all_assignments.append(assignment)

evidence_df = pd.DataFrame([all_genes, all_assignments]).T
evidence_df.columns = ["gene", "evidence"]

# Collect output and write to file
evidence_df = (
    evidence_df.merge(pathway_results_df, how="left", left_on="gene", right_on="symbol")
    .merge(gwas_genes_df, how="left", left_on="gene", right_on="MAPPED_GENE")
    .set_index("gene")
    .drop(["MAPPED_GENE"], axis="columns")
    .sort_index()
    .to_csv(output_file, sep=",", index=True)
)
