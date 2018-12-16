"""
2016 Gregory Way
scripts/summarize_evidence.py

Description:
Takes in genes and evidence support and assigns each gene to the TAD

Usage:
Command line:

     python scripts/summarize_evidence.py

With the following flags:

     --evidence         The location of the evidence file
     --snps             The location of TAD based SNP file
     --group            The group to subset the SNP file (optional)
     --output_file      Where to save the final output

Output:
Trait specific .tsv files of one column each indicating all the genes that fall
in signal TADs
"""

import os
import argparse
import pandas as pd

# Load Command Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-e", "--evidence", help="Location of evidence file")
parser.add_argument("-s", "--snps", help="location of TAD mapped SNPs")
parser.add_argument("-g", "--group", help="Group to subset evidence file",
                    default=None)
parser.add_argument("-o", "--output_file", help="location to write results")
args = parser.parse_args()

# Load Constants
evidence_file = args.evidence
snp_file = args.snps
snp_group = args.group
output_file = args.output_file

# Load data
gene_index = os.path.join('data', 'GENE_index_hg19_hESC.tsv.bz2')
gene_df = pd.read_table(gene_index, index_col=0)

evidence_df = pd.read_csv(evidence_file)

tad_gwas_df = (
    pd.read_csv(snp_file, sep='\t')
    .dropna(subset=['TADidx'])
    .reset_index(drop=True)
    .drop_duplicates('custom_snp')
    )

if snp_group:
    tad_gwas_df = tad_gwas_df.query('group == @snp_group')
    tad_gwas_df = tad_gwas_df.reset_index(drop=True)


def buildTADkey(gwas_snp):
    """
    gwas_snp - a single row subset of the TAD-GWAS input file

    output - The lookup info in the TAD gene dictionary
    i.e. [(Chromosome, TAD_ID:TAD_Start-TAD_End)
    """
    try:
        chrom = gwas_snp['chrom'].replace('chr', '')
    except:
        chrom = gwas_snp['chrom']

    start = int(gwas_snp['TADStart'])
    end = int(gwas_snp['TADEnd'])
    tad_num = int(gwas_snp['TADidx'])
    output = str(tad_num) + ':' + str(start) + '-' + str(end)
    evidence_key = 'chr{}:{}'.format(str(chrom), output)
    return evidence_key


def parse_ev_key(tadkey):
    """
    tadkey - the key for the evidence dictionary

    output - Information stored in the evidence key
    i.e. [ID, Chromosome, Start, End, UCSC lookup]
    """
    chrom, ID, coord = tadkey.split(':')
    start, end = coord.split('-')
    ucsc = '{}:{}-{}'.format(chrom, start, end)
    return [ID, chrom, start, end, ucsc]


# Investigate each significant TAD
evidence_list = []
for tad_row in range(tad_gwas_df.shape[0]):
    # Extract the SNP information
    snp_info = tad_gwas_df.iloc[tad_row, :]

    # Build the key to lookup TAD in dict
    tad_key = buildTADkey(snp_info)

    # Parse the tad key for saving info
    ID, chrom, start, end, ucsc = parse_ev_key(tad_key)

    # Lookup the TAD and extract TAD based genes
    tad_subset_info = gene_df.query('TAD_id == @ID')
    tad_genes = tad_subset_info['gene_name']

    # Subset the evidence dataframe to TAD based genes
    evidence_sub = evidence_df.query('gene in @tad_genes')

    if evidence_sub.shape[0] == 0:
        continue

    # What SNP(s) are associating this TAD?
    snp = '|'.join(set(evidence_sub['snp'].dropna().tolist()))

    # Assign group(s) are we considering
    group = '|'.join(set(evidence_sub['group'].dropna().tolist()))

    # Make sure the SNPs and Groups populate
    evidence_sub = evidence_sub.fillna({'snp': snp, 'group': group})

    # Assign the remaining info to this evidence dataframe
    evidence_sub = evidence_sub.assign(TAD_ID=int(ID),
                                       chromosome=chrom,
                                       TAD_Start=start,
                                       TAD_End=end,
                                       TAD_UCSC=ucsc)
    evidence_list.append(evidence_sub)

# Write out results to file
all_evidence_df = pd.concat(evidence_list, axis='rows')

col = ['gene', 'snp', 'evidence', 'group', 'TAD_ID', 'chromosome', 'TAD_Start',
       'TAD_End', 'TAD_UCSC', 'go_id', 'go_name', 'pval', 'adjP']

(all_evidence_df
    .reindex(col, axis='columns')
    .sort_values(by=['TAD_ID', 'evidence', 'gene'])
    .drop_duplicates()
    .to_csv(output_file, sep='\t', index=False))
