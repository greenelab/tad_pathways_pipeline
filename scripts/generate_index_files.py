"""
2016 Gregory Way
scripts/generate_index_files.py

Description:
Generates index files mapping genomic content to TADs enabling quick lookup

1) SNPs
2) Genes
3) Repetitive Elements

Usage:
Is called by 'scripts/run_pipeline.sh':

      python scripts/generate_index_files.py --TAD-Boundary 'hESC'

Output:
3 gzipped tsv index files and a tsv file of all genes that span TAD boundaries
"""

import os
import argparse
import pandas as pd
from tad_util.util import load_tad, parse_gene_gtf

pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--TAD-Boundary', help='boundary cell type')
parser.add_argument('-f', '--TAD-File', help='path to 3-column Tab-separated TAD domain bed file')
args = parser.parse_args()

TAD_CELL = args.TAD_Boundary

BASE_DIR = 'data/'
GENOME = 'hg19'
REF_GENE_GTF = BASE_DIR + 'gencode.v19.annotation.gtf.gz'
TAD_LOC = args.TAD_File
REF_SNPS = BASE_DIR + 'hg_common-snps.tsv'

REPEAT_FH = BASE_DIR + GENOME + '.fa.out.tsv'
bmd_ld_windows = os.path.join('data', 'BMD_ldwindows.tsv')

SNP_INDEX = 'index/SNP_index_' + GENOME + '_' + TAD_CELL + '.tsv.bz2'
GENE_INDEX = 'index/GENE_index_' + GENOME + '_' + TAD_CELL + '.tsv.bz2'
REPEAT_INDEX = 'index/REPEATS_index_' + GENOME + '_' + TAD_CELL + '.tsv.bz2'
SPANNED_GENES_FH = 'tables/SpannedGenesAcross_' + TAD_CELL + '_TADs.tsv'
bmd_window_genes_file = os.path.join('data', 'BMD_LDwindow_genes.tsv')


def curate_tad_elements(tad_df, input_df, gen_class):
    """
    Loop through TAD boundaries to assign genomic elements to TADs

    Arguments:
    :param tad_df: pandas dataframe, TAD boundary locations
    :param input_df: pandas dataframe, the genomic content to subset
    :param gen_class: str, either 'snp', 'gene', or 'repeat'

    Output:
    Large summary dataframes of the TAD assignment for input genomic element
    """

    big_out_df = pd.DataFrame()
    bnd_df = pd.DataFrame()

    if gen_class in ['gene', 'repeat']:
        input_df['chromosome'] = input_df['chromosome'].map(lambda x: x[3:])

    for tad in tad_df.itertuples():
        tad_id, chrom, start, end = list(tad)
        start, end = int(start), int(end)

        if gen_class == 'snp' and chrom != 'X' and chrom !='Y':
            chrom = int(chrom)

        if gen_class == 'LD':
            chrom = 'chr{}'.format(chrom)

        chm_sub_df = input_df[input_df['chromosome'] == chrom]

        if gen_class == 'snp':
            elem_sub_df = chm_sub_df[(chm_sub_df['position'] >= start) &
                                     (chm_sub_df['position'] < end)]

        else:
            elem_sub_df = chm_sub_df[(chm_sub_df['start'] >= start) &
                                     (chm_sub_df['stop'] <= end)]

            # Get the overlapping genomic elements
            o_df = chm_sub_df[((chm_sub_df['start'] >= start) &
                              (chm_sub_df['start'] < end) &
                              (chm_sub_df['stop'] > end)) |
                              ((chm_sub_df['stop'] >= start) &
                              (chm_sub_df['stop'] < end) &
                              (chm_sub_df['start'] < start))]

            elem_sub_df = pd.concat([elem_sub_df, o_df], axis=0,
                                    ignore_index=True)
            bnd_df = bnd_df.append(o_df, ignore_index=True)

        if gen_class == 'LD':
            # Assign LD information
            elem_sub_df['rs_id'] = tad_id
            elem_sub_df['LD_window_start'] = start
            elem_sub_df['LD_window_end'] = end
        else:
            # Assign TAD information
            elem_sub_df['TAD_id'] = tad_id
            elem_sub_df['TAD_start'] = start
            elem_sub_df['TAD_end'] = end

        big_out_df = big_out_df.append(elem_sub_df, ignore_index=True)

    if gen_class in ['gene', 'repeat']:
        bnd_df = bnd_df.drop_duplicates()

    return big_out_df, bnd_df


# For parsing messy repeat info data
def rm_paren(x): return x.strip('(').strip(')')

# Create index path
if not os.path.exists('index'):
    os.makedirs('index')

# Read in TAD boundary file
tad_df = load_tad(TAD_LOC)

####################################
# PART 1 - SNPs
####################################
snp_df = pd.read_table(REF_SNPS)

big_snp_df, _ = curate_tad_elements(tad_df, snp_df, gen_class='snp')

# Every SNP that did not map is in a boundary region
boundary_snp_df = snp_df.query("rsid not in @big_snp_df.rsid")
boundary_snp_df['TAD_id'] = "Boundary"
boundary_snp_df['TAD_start'] = 0
boundary_snp_df['TAD_end'] = 0

all_snp_df = big_snp_df.append(boundary_snp_df)
all_snp_df.to_csv(SNP_INDEX, sep='\t', compression='bz2')

####################################
# PART 2 - Genes
####################################
# Process data
refGene_df = pd.read_table(REF_GENE_GTF, skiprows=5, header=None)
refGene_df = refGene_df.drop(refGene_df.columns[[5, 7]], axis=1)
refGene_df.columns = ['chromosome', 'db', 'type', 'start', 'stop', 'strand',
                      'info']
parsed_info_df = refGene_df.apply(parse_gene_gtf, axis=1)
parsed_info_df.columns = ['gene_type', 'gene_name']
refGene_df = pd.concat([refGene_df, parsed_info_df], axis=1)
refGene_df = refGene_df[refGene_df['type'] == 'gene'].drop('info', axis=1)

big_gene_df, bound_gene_df = curate_tad_elements(tad_df, refGene_df,
                                                 gen_class='gene')

# Which genes belong entirely in the boundary?
only_bound_gene_df = refGene_df.ix[~refGene_df['gene_name']
                                   .isin(big_gene_df['gene_name'])]
only_bound_gene_df = only_bound_gene_df.ix[~only_bound_gene_df['chromosome']
                                           .isin(['X', 'Y', 'M'])]
only_bound_gene_df['TAD_id'] = "Boundary"
only_bound_gene_df['TAD_start'] = 0
only_bound_gene_df['TAD_end'] = 0

big_gene_tad_df = big_gene_df.append(only_bound_gene_df)
big_gene_tad_df.to_csv(GENE_INDEX, sep='\t', compression='bz2')
bound_gene_df.to_csv(SPANNED_GENES_FH, sep='\t')

####################################
# PART 3 - Repeat Elements #########
####################################
repeats_df = pd.read_csv(REPEAT_FH, delimiter='\t',
                         skiprows=[0, 1, 2], header=None,
                         names=['div', 'chromosome', 'start', 'stop',
                                'repeat'],
                         index_col=False, usecols=[2, 5, 6, 7, 11],
                         converters={'start': rm_paren, 'stop': rm_paren})
repeats_df['start'] = repeats_df['start'].astype(int)
repeats_df['stop'] = repeats_df['stop'].astype(int)

big_rep_tad_df, boundary_rep_df = curate_tad_elements(tad_df, repeats_df,
                                                      gen_class='repeat')

big_rep_tad_df = big_rep_tad_df.append(boundary_rep_df, ignore_index=True)
big_rep_tad_df.to_csv(REPEAT_INDEX, sep='\t', compression='bz2')

####################################
# PART 4 - BMD genes in LD Windows #
####################################
bmd_ld_df = pd.read_table(bmd_ld_windows, index_col=0)
bmd_ld_genes_df, _ = curate_tad_elements(bmd_ld_df, refGene_df, gen_class='LD')
bmd_ld_genes_df.to_csv(bmd_window_genes_file, sep='\t', index_col=0)
