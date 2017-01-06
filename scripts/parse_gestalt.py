"""
2016 Gregory Way
scripts/parse_gestalt.py

Description:
Take as input the .tsv results from the WebGestalt analysis and outputs trait
specific TAD pathway genes

Usage:
Command line 'python scripts/parse_gestalt.py' with the required flag:

    --trait     is predefined by the manual WebGestalt step (see README). Or it
                can be specified by the name of the gestalt output file name.

                    Example 1: "bmd" for Bone Mineral Density GWAS
                    Example 2: "custom_group" for file named -
                               "gestalt/custom_group_gestalt.tsv"

And the optional flag:

    --process_type      include flag if using previous version of WebGestalt

Output:
csv file with collapsed pathway specific information
csv file with summary of pathway p values
"""

import os
import argparse
import pandas as pd
import io
import re


def read_gestalt(file_name, process):
    """
    Read in the gestalt data and parse pathways

    Arguments:
    :param file_name: the location of the file to read
    :param process: the type of Gestalt file to parse

    Output:
    1) Pandas Dataframe of entire TRAIT pathways
    2) Pathway summary text file
    """

    # Read in the entire file as string
    with open(file_name, 'r') as trait_fh:
        gestalt = trait_fh.read()

    if process:
        gestalt = re.sub('\n\t\t\t\t\t\n', '\n\n', gestalt)
        gestalt = re.sub('\t\t\t\t\t\n', '\n', gestalt)
        gestalt = re.sub('\t\t\t\n', '\n', gestalt)

    # Remove trailing whitespace from file
    gestalt = gestalt.rstrip()

    # Triple new lines demarcate pathways, skip header info
    stanzas = gestalt.split('\n\n\n')[1:]

    # Initialize empty list to append pathway DataFrames
    stanzas_df_list = list()

    # Loop through each pathway and combine into large DataFrame
    pathway_names = []
    pathway_adjp = []
    pathway_numgenes = []
    for pathway in stanzas:
        pathway_file = io.StringIO(pathway)

        # Parse pathway and enrichment information
        go_info = next(pathway_file).rstrip('\n').split('\t')
        enrichment = next(pathway_file).rstrip('\n').split(';')
        go_class, go_name, go_id = go_info
        adj_p = float(enrichment[5][5:])

        # Read in the io file
        pathway_header = ['symbol', 'name', 'entrez_id', 'ensemble_id']
        pathway_df = pd.read_table(pathway_file, skip_blank_lines=True,
                                   names=pathway_header,
                                   usecols=[0, 3, 4, 5])

        # Add GO information to the pathway DataFrame
        pathway_df['go_class'] = go_class
        pathway_df['go_name'] = go_name
        pathway_df['go_id'] = go_id
        pathway_df['adjP'] = adj_p
        stanzas_df_list.append(pathway_df)

        # Add to general information
        pathway_names.append(go_name)
        pathway_adjp.append(adj_p)
        pathway_numgenes.append(pathway_df.shape[0])

    return_df = pd.concat(stanzas_df_list)
    return_df = return_df.sort_values(by='adjP')
    gen_info = return_df.groupby(['go_id', 'go_name',
                                  'adjP']).symbol.count().reset_index()
    gen_info.rename(columns={'symbol': 'num_genes'}, inplace=True)
    gen_info = gen_info.sort_values(by=['adjP', 'go_name'])

    return return_df, gen_info

if __name__ == '__main__':

    # Load Command Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--trait', help='symbol for trait data')
    parser.add_argument('-p', '--process_type',
                        help='deals with different Gestalt file types',
                        action='store_false')
    args = parser.parse_args()

    # Load Constants
    trait = args.trait
    process = args.process_type

    trait_file = os.path.join('gestalt', '{}_gestalt.tsv'.format(trait))
    trait_out_file = os.path.join('gestalt',
                                  '{}_complete_gestalt.tsv'.format(trait))
    out_p_val = os.path.join('gestalt', '{}_pvalues.tsv'.format(trait))

    # Load Data
    gestalt_data, gestalt_p = read_gestalt(trait_file, process)

    # Write files
    gestalt_data.to_csv(trait_out_file, sep='\t', index=False)
    gestalt_p.to_csv(out_p_val, sep='\t', index=False)
