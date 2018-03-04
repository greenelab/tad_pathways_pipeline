"""
2016 Gregory Way
scripts/summarize_evidence.py

Description:
Takes in genes and evidence support and assigns each gene to the TAD

Usage:
Command line:

     python scripts/assign_evidence_to_TADs.py

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
import csv

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

tad_gwas_df = pd.read_csv(snp_file, sep='\t')
tad_gwas_df = tad_gwas_df.dropna(subset=['TADidx'])
tad_gwas_df = tad_gwas_df.reset_index(drop=True)
if snp_group:
    tad_gwas_df = tad_gwas_df[tad_gwas_df['group'] == snp_group]
    tad_gwas_df = tad_gwas_df.reset_index(drop=True)


def buildTADkey(gwas_snp):
    """
    gwas_snp - a single row subset of the TAD-GWAS input file

    output - The lookup info in the TAD gene dictionary
    i.e. [(Chromosome, TAD_ID:TAD_Start-TAD_End)
    """

    chrom = gwas_snp['chrom'].replace('chr', '') 

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
    chrom = tadkey.split(':')
    ID = chrom[1]
    start = chrom[2].split('-')
    end = start[1]
    start = start[0]
    chrom = chrom[0]
    ucsc = chrom + ':' + start + '-' + end
    return [ID, chrom, start, end, ucsc]

# Investigate each significant TADs
evidence_dict = {}
for tad_row in range(len(tad_gwas_df)):
    snp_info = tad_gwas_df.ix[tad_row, :]

    # Build the key to lookup TAD in dict and lookup
    e_key = buildTADkey(snp_info)

    if e_key not in evidence_dict.keys():
        evidence_dict[e_key] = []

    tad_index = int(snp_info['TADidx'])
    tad_subset_info = gene_df.ix[gene_df['TAD_id'] == str(tad_index), :]
    tad_genes = tad_subset_info['gene_name'].tolist()

    # Subset the evidence dataframe to TAD based genes
    evidence_sub = evidence_df.ix[evidence_df['gene'].isin(tad_genes), :]

    # Loop over each of the rows
    TAD_gene_evidence = []
    for ev in range(0, evidence_sub.shape[0]):
        gene = evidence_sub['gene'].tolist()[ev]
        ev_type = evidence_sub['evidence'].tolist()[ev]
        TAD_gene_evidence.append([gene, ev_type])

    evidence_dict[e_key].append(TAD_gene_evidence)

# Write out results to file
with open(output_file, 'w') as out_fh:
    tadwriter = csv.writer(out_fh, delimiter='\t')
    tadwriter.writerow(['gene', 'evidence', 'TAD ID', 'chromosome',
                        'TAD Start', 'TAD End', 'UCSC'])
    for tadkey in evidence_dict.keys():
        ID, chrom, start, end, ucsc = parse_ev_key(tadkey)
        for gene_evidence_list in evidence_dict[tadkey][0]:
            gene, evidence = gene_evidence_list
            tadwriter.writerow([gene, evidence, ID, chrom, start, end, ucsc])
