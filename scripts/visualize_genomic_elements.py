"""
2016 Gregory Way
scripts/visualize_genomic_elements.py

Description:
Summarizes the location of genomic elements across TADs

Usage:
Is called by 'scripts/visualize.sh' which is run inside of
'scripts/run_pipeline.sh'. This particular script will output the location
of genomic elements in a given input TAD

      python scripts/visualize_genomic_elements.py --TAD-Boundary 'hESC'

Output:
Several .pdf plots in "figures/genome/" and chisquare analyses of the
"rightness" of SNPs in TADs and protein coding genes near boundaries.
"""

import os
import argparse

import csv
import numpy as np
import pandas as pd
from scipy.stats import chisquare
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

from tad_util.util import assign_bin

plt.figure.max_open_warning = 0
sns.set_style("whitegrid")
sns.set_style("ticks")
sns.set_context("paper", rc={"font.size": 20, "axes.titlesize": 20,
                             "axes.labelsize": 20, "xtick.labelsize": 12,
                             "ytick.labelsize": 12})

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--TAD-Boundary', help='boundary cell type') 
args = parser.parse_args()

# Load Constants
num_bins = 50
tad_cell = args.TAD_Boundary
xlab = [''] * num_bins
for x in range(0, 50, 10):
    xlab[x] = x

genome = 'hg19'

# Input files
base_file = '{}_{}'.format(genome, tad_cell)

snp_index = os.path.join('index', 'SNP_index_{}.tsv.bz2'.format(base_file))
gene_index = os.path.join('index', 'GENE_index_{}.tsv.bz2'.format(base_file))
repeat_index = os.path.join('index', 'REPEATS_index_{}.tsv.bz2'
                                     .format(base_file))

# Output files
fig_base = 'figures'

if not os.path.exists(fig_base):
    os.makedirs(fig_base)

snp_count_file = os.path.join(fig_base, 'snp_count_{}.pdf'.format(base_file))
snp_dist_file = os.path.join(fig_base, 'snp_tad_distribution_{}.pdf'
                                       .format(base_file))
snp_chrom_file = os.path.join(fig_base, 'snp_tad_distrib_chromosomes_{}.pdf'
                                        .format(base_file))
snp_chi_square = os.path.join('results',
                              'tad_snp_rightness_chi_{}.csv').format(base_file)

gene_count_file = os.path.join(fig_base, 'gene_count_{}.pdf'
                                         .format(base_file))
gene_chrom_file = os.path.join(fig_base, 'gene_tad_distrib_chromosomes_{}.pdf'
                                         .format(base_file))
gene_type_file = os.path.join(fig_base, 'gene_types_{}.pdf'.format(base_file))
gene_chi_square = os.path.join('results',
                               'tad_gene_bound_chi_{}.csv').format(base_file)

repeat_count_file = os.path.join(fig_base, 'repeat_count_{}.pdf'
                                           .format(base_file))
rep_type_file = os.path.join(fig_base, 'repeat_type_{}_.pdf'.format(base_file))
repeat_dist = os.path.join(fig_base, 'repeat_type_all_distrib_{}.pdf'
                                     .format(base_file))

# Load Data
gene_types_df = pd.read_table(os.path.join('tables',
                                           'gene_classification.tsv'))
snp_df = pd.read_table(snp_index, index_col=0)
gene_df = pd.read_table(gene_index, index_col=0)
repeat_df = pd.read_table(repeat_index, index_col=0)

#########################
# PART 1 - SNPs
#########################
# Process SNP dataframe
snp_df = snp_df[snp_df['TAD_id'] != 'Boundary']
bin_s = snp_df.apply(lambda x: assign_bin(x, bins=num_bins, ID='SNP'), axis=1)
snp_df = snp_df.assign(tad_bin=bin_s)

# Jointplot of number of SNPs per TAD by TAD length
plot_ready = snp_df.assign(tad_length=np.log10(snp_df.TAD_end
                                                     .sub(snp_df.TAD_start)))
plot_ready = pd.DataFrame(plot_ready.groupby(['TAD_id', 'tad_length'])
                                    .tad_bin.count()).reset_index()
plot_ready = plot_ready.assign(snp_count_alt=plot_ready.tad_bin.div(1000))
ax = sns.jointplot('tad_length', 'snp_count_alt', data=plot_ready,
                   kind='scatter', stat_func=None,
                   color=sns.xkcd_rgb['medium green'], joint_kws={'s': 3})
ax.set_axis_labels(xlabel='TAD Length (log10 kb)',
                   ylabel='Number of SNPs (x1000)')
plt.tight_layout()
plt.savefig(snp_count_file)
plt.close()

# Distribution of SNPs across TADs
summary_snp = snp_df['tad_bin'].value_counts(sort=False)
p = sns.pointplot(x=summary_snp.index, y=summary_snp / 1000,
                  color=sns.xkcd_rgb["medium green"], scale=0.5)
sns.despine()
p.set(xticklabels=xlab)
p.set(ylabel='Number of SNPs (x1000)', xlabel='TAD Bins')
p.set_title('Distribution of SNPs across TADs')
plt.tight_layout()
plt.savefig(snp_dist_file)
plt.close()

# Chromosome-specific distribution
snp_chrom = snp_df.groupby('chromosome').tad_bin.value_counts(sort=False).\
  unstack(level=0)
with PdfPages(snp_chrom_file) as pdf:
    for chrom, chrom_df in snp_chrom.iteritems():
        p = sns.pointplot(x=chrom_df.index, y=chrom_df,
                          color=sns.xkcd_rgb["medium green"], scale=0.5)
        sns.despine()
        p.set(xticklabels=xlab)
        p.set(ylabel='Number of SNPs', xlabel='TAD Bins')
        p.set_title('SNP Distribution in Chromosome {}'.format(chrom))
        plt.tight_layout()
        pdf.savefig()
        plt.close()

# SNPs appear to be more concentrated on the right side of TADs
snp_side = [snp_df[snp_df['tad_bin'] < 25].shape[0],
            snp_df[snp_df['tad_bin'] >= 25].shape[0]]
tad_snp_sig = chisquare(snp_side)

with open(snp_chi_square, 'w') as chisq_fh:
    snpwriter = csv.writer(chisq_fh, delimiter=',')
    snpwriter.writerow(['SNPs in the left vs. right of {} TAD'
                        .format(tad_cell)])
    snpwriter.writerow(['left', 'right'])
    snpwriter.writerow(snp_side)
    snpwriter.writerow(tad_snp_sig)

#########################
# PART 2 - Genes
#########################
# Process genes
gene_df = gene_df[gene_df['TAD_id'] != 'Boundary']
bin_assign_gene = gene_df.apply(lambda x: assign_bin(x, bins=num_bins,
                                                     ID='gene'), axis=1)
gene_df = gene_df.assign(tad_bin=bin_assign_gene)
gene_df = gene_df[gene_df['tad_bin'] != -1]

# Jointplot of number of Genes per TAD
plot_ready_gene = gene_df.assign(tad_length=np.log10(gene_df.TAD_end
                                                     .sub(gene_df.TAD_start)))
plot_ready_gene = pd.DataFrame(plot_ready_gene.groupby(['TAD_id',
                                                        'tad_length'])
                                              .tad_bin.count()).reset_index()
plot_ready_gene = plot_ready_gene.assign(gene_count_alt=plot_ready_gene
                                         .tad_bin)
ax = sns.jointplot('tad_length', 'gene_count_alt', data=plot_ready_gene,
                   kind='scatter', stat_func=None,
                   color=sns.xkcd_rgb['medium green'], joint_kws={'s': 3})
ax.set_axis_labels(xlabel='TAD Length (log10 kb)',
                   ylabel='Number of Genes')
plt.savefig(gene_count_file)
plt.close()

# Chromosome specific distribution of genes across TADs
gene_chrom = gene_df.groupby('chromosome').tad_bin.value_counts(sort=False).\
  unstack(level=0)
with PdfPages(gene_chrom_file) as pdf:
    for chrom, chrom_df in gene_chrom.iteritems():
        ax = sns.pointplot(x=chrom_df.index, y=chrom_df,
                           color=sns.xkcd_rgb["medium green"], scale=0.5)
        sns.despine()
        ax.set(xticklabels=xlab)
        ax.set(ylabel='Number of Genes', xlabel='TAD Bins')
        ax.set_title('Gene Distribution in Chromosome {}'.format(chrom))
        plt.tight_layout()
        pdf.savefig()
        plt.close()

# Gene-type specific distribution across TADs
gene_types_df = gene_types_df[gene_types_df[genome] == 1]
summary_gene_classes = []
with PdfPages(gene_type_file) as pdf:
    for idx, gene in gene_types_df.iterrows():
        gene_class = gene['gene_class']
        gene_type = gene['gene_type']

        if gene_class in ['tr_gene', 'ig_gene', 'tr_pseud', 'ig_pseud']:
            gene_type = gene_types_df[gene_types_df['gene_class'] ==
                                      gene_class]['gene_type']
            gene_sub_df = gene_df[gene_df['gene_type'].isin(gene_type)]
            plot_title = gene_class

            if gene_class in summary_gene_classes:
                continue
            else:
                summary_gene_classes.append(gene_class)

        elif gene_class == 'std' and gene_type != 'all':
            gene_sub_df = gene_df[gene_df['gene_type'] == gene_type]
            plot_title = gene_type
        elif gene_type == 'all':
            gene_sub_df = gene_df
            plot_title = 'Distribution of Genes across TADs'

        sum_gene = gene_sub_df['tad_bin'].value_counts(sort=False).sort_index()

        ax = sns.pointplot(x=sum_gene.index, y=sum_gene,
                           color=sns.xkcd_rgb["medium green"], scale=0.5)
        sns.despine()
        ax.set(xticklabels=xlab)
        ax.set(ylabel='Number of Genes', xlabel='TAD Bins')
        ax.set_title(plot_title)
        plt.tight_layout()
        pdf.savefig()
        plt.close()

# Chisquare of genes on TAD boundaries
protein_coding = gene_df[gene_df['gene_type'] == 'protein_coding']
bin_list = list(range(num_bins))[0:2] + list(range(num_bins))[-2:]
boundary_df = protein_coding[protein_coding['tad_bin'].isin(bin_list)]
num_genes_b = boundary_df.shape[0]
num_genes_c = protein_coding.shape[0] - num_genes_b
chi_test = [num_genes_b, num_genes_c]
exp = protein_coding.shape[0] / num_bins
bound_chi = chisquare(chi_test, f_exp=[exp * len(bin_list),
                                       exp * (num_bins - len(bin_list))])

with open(gene_chi_square, 'w') as chisq_fh:
    genewriter = csv.writer(chisq_fh, delimiter=',')
    genewriter.writerow(['Genes at boundaries vs. center of {} TAD'
                        .format(tad_cell)])
    genewriter.writerow(['bound', 'center'])
    genewriter.writerow(chi_test)
    genewriter.writerow(bound_chi)

#########################
# PART 3 - Repeats
#########################
# Process Repeats
repeat_df = repeat_df.fillna('Boundary')
repeat_df = repeat_df[repeat_df['TAD_id'] != 'Boundary']
bin_assign_repeat = repeat_df.apply(lambda x: assign_bin(x, bins=num_bins,
                                                         ID='repeat'), axis=1)
repeat_df = repeat_df.assign(tad_bin=bin_assign_repeat)
repeat_df = repeat_df[repeat_df['tad_bin'] != -1]

# Jointplot of number of repeats per TAD
repeat_df.TAD_end = repeat_df.TAD_end.astype(int)
repeat_df.TAD_start = repeat_df.TAD_start.astype(int)
plot_ready_repeat = repeat_df.assign(tad_length=np.log10(repeat_df.TAD_end
                                     .sub(repeat_df.TAD_start)))
plot_ready_repeat = pd.DataFrame(plot_ready_repeat.groupby(['TAD_id',
                                                            'tad_length'])
                                 .tad_bin.count()).reset_index()
plot_ready_repeat = plot_ready_repeat.assign(rep_count_alt=plot_ready_repeat
                                             .tad_bin.div(100))
ax = sns.jointplot('tad_length', 'rep_count_alt', data=plot_ready_repeat,
                   kind='scatter', stat_func=None,
                   color=sns.xkcd_rgb['medium green'], joint_kws={'s': 3})
ax.set_axis_labels(xlabel='TAD Length (log10 kb)',
                   ylabel='Number of Repeats (x100)')
plt.savefig(repeat_count_file)
plt.close()

# Distribution of different classes of repeats across TADs
with PdfPages(rep_type_file) as pdf:
    for repeat_type in repeat_df['repeat'].unique():
        if '?' not in repeat_type:
            repeat_fh = repeat_type.replace('/', '_')
            rep_sub = repeat_df[repeat_df['repeat'] == repeat_type]
            sum_rep = rep_sub['tad_bin'].value_counts(sort=False).sort_index()
            p = sns.pointplot(x=sum_rep.index, y=sum_rep,
                              color=sns.xkcd_rgb["medium green"], scale=0.5)
            sns.despine()
            p.set(xticklabels=xlab)
            p.set(ylabel='Number of Repeats', xlabel='TAD Bins')
            p.set_title(repeat_type + ' Distribution')
            plt.tight_layout()
            pdf.savefig()
            plt.close()

# Distribution of all repeats
sum_repeat = repeat_df['tad_bin'].value_counts(sort=False).sort_index()
p = sns.pointplot(x=sum_repeat.index, y=sum_repeat.div(100),
                  color=sns.xkcd_rgb["medium green"], scale=0.5)
sns.despine()
p.set(xticklabels=xlab)
p.set(ylabel='Number of Repeats (x100)', xlabel='TAD Bins')
p.set_title('All Repeats Distribution')
plt.tight_layout()
plt.savefig(repeat_dist)
plt.close()
