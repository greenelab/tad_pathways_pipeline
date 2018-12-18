"""
(C) 2016 Gregory Way
visualize_TAD_locations.py

Description:
Observe SNP/gene/repeat distributions and locations across TADs (hg19)

Usage:
Is called by 'ANALYSIS.sh'.

    python visualize_TAD_locations.py --genome 'hg19' --TAD-Boundary 'hESC'

Where flags for --TAD-Boundary are: "hESC", "IMR90", "mESC", or "cortex"
and the flags for --genome are: "hg19" or "mm9"

Output:
1) SNP distribution and histogram across TADs
2) Gene distribution
    A) Several matplotlib plots in the 'figures/' folder
    B) Several gene location plots across TADs in 'figures/gene/'
    C) Chisquare analysis results for gene start sites near boundaries
3) Repeat distribution
    A) Several matplotlib plots in the 'figures/' folder
    B) Several repeat location plots across TADs in 'figures/repeats/'
"""

import argparse
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
from scipy.stats import chisquare

from tad_util.util import parse_TAD_name, parse_SNP_position, parse_gene_info
from tad_util.util import parse_repeat_info, ID_TAD_bins

plt.figure.max_open_warning = 0

# Load Command Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--genome', help='Can be either hg19 or mm9')
parser.add_argument('-t', '--TAD-Boundary', help='boundary cell type. The'
                    'options can be "hESC", "IMR90", "mESC", or "cortex"')
args = parser.parse_args()

####################################
# Load Constants
####################################
NUM_BINS = 50
TAD_CELL = args.TAD_Boundary
GENOME = args.genome

# Generate file names depending on input
if TAD_CELL in ['hESC', 'IMR90', 'mESC', 'cortex']:
    if GENOME in ['hg19', 'mm9']:
        # Files to load
        INDEX_BASE = GENOME + '_' + TAD_CELL + '.p'
        SNP_INDEX = 'index/SNP_index_' + INDEX_BASE
        GENE_INDEX = 'index/GENE_index_' + INDEX_BASE
        REPEAT_INDEX = 'index/REPEAT_index_' + INDEX_BASE

        # Consistent file name conventions
        FIG_BASE = 'figures/' + GENOME + '/'
        FIG_END = GENOME + '_' + TAD_CELL + '.png'

        # Files to save
        TAD_SNP_CHI = 'results/TAD_SNP_rightness_Chisquare_' + \
                      GENOME + '_' + TAD_CELL + '.txt'
        SNP_DIST = FIG_BASE + 'SNPdist_' + FIG_END
        SNP_CHROM = FIG_BASE + 'chrom/SNPdistribution_chrom_'
        SNP_HIST = FIG_BASE + 'SNPcounts_' + FIG_END
        GENE_TAD_BOUND = 'tables/genes_at_' + TAD_CELL + '_TAD_boundaries.tsv'
        GENE_CHI = 'results/' + TAD_CELL + '_TADBoundary_Chisquare.txt'
        GENE_LINE = FIG_BASE + 'gene/genelocation_'
        GENE_CHROM = FIG_BASE + 'chrom/GeneDistribution_chrom_'
        GENE_ALL_CHROM = FIG_BASE + 'chrom/GeneDist_chromo' + FIG_END
        GENE_HIST = FIG_BASE + TAD_CELL + '_TAD_Gene_histogram.png'
        TAD_LEN_HIST = FIG_BASE + TAD_CELL + '_TAD_length_histogram.png'
        REP_TYPE = FIG_BASE + 'repeats/repeat_' + TAD_CELL + '_'
        REP_CHROM = FIG_BASE + 'chrom/RepeatDist_chrom_'
        REP_ALL_CHROM = FIG_BASE + 'chrom/RepeaDist_chrom_' + FIG_END
        REP_HIST = FIG_BASE + TAD_CELL + '_TAD_Repeat_histogram.png'
    else:
        raise ValueError('Please input either "hg19" or "mm9"')

else:
    raise ValueError('Please input either "hESC", "IMR90", "mESC", or '
                     '"cortex"')

# Gene types described http://www.gencodegenes.org/gencode_biotypes.html
gene_types = ['protein_coding', 'pseudogene', 'miRNA', 'snRNA', 'snoRNA',
              'lincRNA', 'rRNA', 'all', 'misc_RNA', 'sense_intronic',
              'processed_transcript', 'sense_overlapping', 'IG_pseudogene',
              'TR_pseudogene', 'TR_gene', 'polymorphic_pseudogene',
              'antisense', 'IG_gene']

TR_genes = ['TR_V_gene', 'TR_J_gene', 'TR_D_gene', 'TR_C_gene']
IG_genes = ['IG_V_gene', 'IG_C_gene', 'IG_C_gene', 'IG_J_gene',
            'IG_D_gene']
TR_pseud = ['TR_V_pseudogene', 'TR_J_pseudogene']
IG_pseud = ['IG_V_pseudogene', 'IG_C_pseudogene', 'IG_J_pseudogene']

if GENOME == 'mm9':
    gene_types = gene_types + ['processed_pseudogene', 'TEC',
                               'transcribed_processed_pseudogene',
                               'unprocessed_pseudogene', 'scaRNA',
                               'transcribed_unprocessed_pseudogene',
                               'ribozyme', 'unitary_pseudogene',
                               '3prime_overlapping_ncRNA',
                               'bidirectional_promoter_lncRNA', 'scRNA',
                               'translated_unprocessed_pseudogene',
                               'transcribed_unitary_pseudogene', 'sRNA',
                               'IG_LV_gene', 'IG_D_pseudogene', 'macro_lncRNA']
else:
    # Gene types described http://www.gencodegenes.org/gencode_biotypes.html
    gene_types = gene_types + ['3prime_overlapping_ncrna']

####################################
# Load Data
####################################
KG_TADs = pd.read_pickle(SNP_INDEX)
TADdictgenes = pd.read_pickle(GENE_INDEX)
TADdictrepeats = pd.read_pickle(REPEAT_INDEX)

####################################
# Set matplotlib defaults
####################################
fig_size = 4, 3
margin = 2.5
plt.rcParams['figure.figsize'] = fig_size
x0, x1, y0, y1 = plt.axis()
fig, ax = plt.subplots()
#plt.axis((x0 + margin, x1 + margin, y0 + margin, y1 + margin))
plt.tight_layout()

####################################
# PART 1 - SNPs ####################
####################################

####################################
# Observe location of SNPs in TADs
####################################
# Initialize a dictionary that will store location information for gene types
TADdict = range(NUM_BINS)
TADdict = dict.fromkeys(TADdict)
for key in TADdict.iterkeys():
    TADdict[key] = 0

# Initialize a chromosome specific dictionary holding bin specific info
Chromosome_Specific_TAD = dict.fromkeys(KG_TADs.keys())
for chromkey in Chromosome_Specific_TAD.keys():
    # Create a list of zeros of len(NUM_BINS)
    Chromosome_Specific_TAD[chromkey] = [0] * (NUM_BINS)

# Also observe distribution of amount of SNPs in each TAD
howmanySNPsinTAD = []

# Loop through each TAD and SNP to fill the histogram
SNP_side = [0, 0]
for key in KG_TADs.iterkeys():
    print key
    if key in ['NoTAD', 'Boundary', 'X', 'Y']:
        continue
    for tadkey, tadval in KG_TADs[key].iteritems():
        howmanySNPsinTAD.append(len(tadval))
        if len(tadval) > 0:
            tadname = parse_TAD_name(tadkey)
            for snp in range(len(tadval)):
                snp_loc = parse_SNP_position(tadval.ix[snp, :])
                bin_assign = ID_TAD_bins(tadname, NUM_BINS, snp_loc)
                TADdict[bin_assign] += 1

                # Count the number of SNPs in left/right
                if bin_assign < 25:
                    SNP_side[0] += 1
                else:
                    SNP_side[1] += 1

                # Add a count for chromsome specificity
                Chromosome_Specific_TAD[key][bin_assign] += 1

# Convert to DataFrame
SNPLocations = pd.DataFrame.from_dict(TADdict, orient='index')
SNPLoc_chrom = pd.DataFrame.from_dict(Chromosome_Specific_TAD, orient='index')
SNPLoc_chrom = SNPLoc_chrom.transpose()

####################
# Perform ChiSquare test of SNPs to 'right' of TAD
####################
TAD_SNP_sig = chisquare(SNP_side)

with open(TAD_SNP_CHI, 'w') as chisq_fh:
    chisq_fh.write('SNPs in the left vs. right of ' + TAD_CELL + ' TAD\n')
    chisq_fh.write('left,right\n')
    chisq_fh.write(','.join('%s' % x for x in SNP_side) + '\n')
    chisq_fh.write(','.join('%s' % x for x in TAD_SNP_sig))

####################################
# Visualization
####################################
# SNP Distribution
if GENOME == 'mm9':
    fig_size = 4, 3
    margin = 10
    plt.rcParams['figure.figsize'] = fig_size
    x0, x1, y0, y1 = plt.axis()
    #plt.axis((x0 + margin, x1 + margin, y0 + margin, y1 + margin))
    font_size = 8
    plt.rcParams['font.size'] = font_size
else:
    font_size = 10
    plt.rcParams['font.size'] = font_size

ax = SNPLocations.plot(legend=False)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
plt.xlabel('Bins (Normalized TAD length)')
plt.ylabel('Frequency')
plt.title('SNP location inside ' + TAD_CELL + ' TADs')
plt.grid(False)
plt.tick_params(axis='both', which='both', bottom='off', left='off', top='off',
                right='off')
plt.savefig(SNP_DIST)
plt.close()

#####################################
# Plot Chromosome Specific Line Graphs
#####################################
for chrom in Chromosome_Specific_TAD.keys():
    if chrom in ['NoTAD', 'Boundary', 'X', 'Y']:
        continue
    nameplot = SNP_CHROM + chrom + '_' + GENOME + '_' + TAD_CELL + '.png'
    chrom_loc = SNPLoc_chrom[chrom]
    chrom_loc.plot()
    plt.xlabel('Bins (Normalized TAD length)')
    plt.ylabel('Frequency')
    plt.title('SNP location inside ' + TAD_CELL + ' TADs\nChromosome ' +
              str(chrom))
    plt.grid(False)
    plt.savefig(nameplot)
    plt.close()

# Histogram
if GENOME == 'mm9':
    margin = 10
    x0, x1, y0, y1 = plt.axis()
    #plt.axis((x0 + margin, x1 + margin, y0 + margin, y1 + margin))
    font_size = 8
    plt.rcParams['font.size'] = font_size
else:
    plt.rcParams['font.size'] = 9

n, bins, patches = plt.hist(howmanySNPsinTAD, NUM_BINS, normed=1,
                            facecolor='green', alpha=0.75)
plt.xlabel('Number of SNPs')
plt.ylabel('Frequency')
plt.title('Distribution of number of SNPs in ' + TAD_CELL + ' TADs')
addtext = 'mean = ' + str(round(np.mean(howmanySNPsinTAD), 2)) + '\n' + \
          'std = ' + str(round(np.std(howmanySNPsinTAD), 2))
if GENOME == 'mm9':
    plt.axis([min(howmanySNPsinTAD), max(howmanySNPsinTAD), 0, 0.0001])
    plt.text(max(howmanySNPsinTAD), .00004, addtext, verticalalignment='top',
             horizontalalignment='right')
else:
    plt.axis([min(howmanySNPsinTAD), max(howmanySNPsinTAD), 0, 0.0005])
    plt.text(max(howmanySNPsinTAD), .0003, addtext, verticalalignment='top',
             horizontalalignment='right')
plt.grid(False)
plt.savefig(SNP_HIST)
plt.close()

####################################
# PART 2 - Genes ###################
####################################
if GENOME == 'mm9':
    fig_size = 4, 3
    margin = 2.5
    plt.rcParams['figure.figsize'] = fig_size
    x0, x1, y0, y1 = plt.axis()
    #plt.axis((x0 + margin, x1 + margin, y0 + margin, y1 + margin))

####################################
# Observe location of genes in TADs
####################################
# Initialize a dictionary that will store location information for gene types
TADgenepos = dict.fromkeys(gene_types)
for genekey in gene_types:
    # Create a list of zeros of len(NUM_BINS)
    TADgenepos[genekey] = [0] * (NUM_BINS)

# Initialize a chromosome specific dictionary holding bin specific info
Chromosome_Specific_TAD = dict.fromkeys(TADdictgenes.keys())
for chromkey in Chromosome_Specific_TAD.keys():
    if chromkey in ['Boundary']:
        continue
    # Create a list of zeros of len(NUM_BINS)
    Chromosome_Specific_TAD[chromkey] = [0] * (NUM_BINS)

# Determine significance of 'leftness' of TADs
TAD_near_boundary = [0, 0]
TAD_nb_protein = [0, 0]
TAD_nb_rRNA = [0, 0]
left_bin = NUM_BINS / 2

# Also, we're interested in the genes near TAD boundaries
gene_TAD_boundary = []
# Obtain the gene position information
for key, value in TADdictgenes.iteritems():
    if key in ['Y', 'Boundary']:
        continue
    for tadkey, tadval in TADdictgenes[key].iteritems():
        tadkey = parse_TAD_name(tadkey)
        for gene in range(len(tadval)):
            # gene_info stores [gene_type, chromosome, start_pos]
            gene_info = parse_gene_info(tadval.ix[gene, :])
            gene_name = tadval.ix[gene, :]['gene']

            binID = ID_TAD_bins(tadkey, NUM_BINS, gene_info[2],
                                IDtype='gene')
            if binID == -1:
                continue

            # Add the bin number to the appropriate dictionary
            gene_type = gene_info[0]
            if gene_type in TR_genes:
                TADgenepos['TR_gene'][binID] += 1
            elif gene_type in IG_genes:
                TADgenepos['IG_gene'][binID] += 1
            elif gene_type in TR_pseud:
                TADgenepos['TR_pseudogene'][binID] += 1
            elif gene_type in IG_pseud:
                TADgenepos['IG_pseudogene'][binID] += 1
            else:
                TADgenepos[gene_type][binID] += 1

            # Add the bin number of all genes
            TADgenepos['all'][binID] += 1

            # Determine which half of bin gene belongs in and store
            if binID in [0, 1, 48, 49]:
                if GENOME == 'hg19':
                    if gene_type == 'protein_coding':
                        TAD_nb_protein[0] += 1
                    elif gene_type == 'rRNA':
                        TAD_nb_rRNA[0] += 1
                TAD_near_boundary[0] += 1
                gene_TAD_boundary.append([gene_name, binID])

            else:
                if GENOME == 'hg19':
                    if gene_type == 'protein_coding':
                        TAD_nb_protein[1] += 1
                    elif gene_type == 'rRNA':
                        TAD_nb_rRNA[1] += 1
                TAD_near_boundary[1] += 1

            # Add a count for chromsome specificity
            Chromosome_Specific_TAD[key][binID] += 1

# Convert to pandas dataframe
geneLocations = pd.DataFrame(TADgenepos)
chromLocations = pd.DataFrame(Chromosome_Specific_TAD)
gene_TAD_boundary = pd.DataFrame(gene_TAD_boundary, columns=['gene', 'bin'])

# Write gene TAD boundary to csv
gene_TAD_boundary.to_csv(GENE_TAD_BOUND, sep='\t')

####################
# Perform ChiSquare test of TAD near boundary
####################
# Expected values for each bin
ex_all = sum(TAD_near_boundary) / NUM_BINS
TAD_nb_sig = chisquare(TAD_near_boundary, f_exp=[ex_all * 4, ex_all * 46])
if GENOME == 'hg19':
    ex_pro = sum(TAD_nb_protein) / NUM_BINS
    ex_rRNA = sum(TAD_nb_rRNA) / NUM_BINS
    TAD_nb_sig_pro = chisquare(TAD_nb_protein, f_exp=[ex_pro * 4, ex_pro * 46])
    TAD_nb_sig_rRNA = chisquare(TAD_nb_rRNA, f_exp=[ex_rRNA * 4, ex_rRNA * 46])

with open(GENE_CHI, 'w') as chisq_fh:
    chisq_fh.write(TAD_CELL + ' near boundaries of all genes\n')
    chisq_fh.write(','.join('%s' % x for x in TAD_near_boundary) + '\n')
    chisq_fh.write(','.join('%s' % x for x in TAD_nb_sig) + '\n')
    if GENOME == 'hg19':
        chisq_fh.write(TAD_CELL + ' near boundaries of protein coding genes\n')
        chisq_fh.write(','.join('%s' % x for x in TAD_nb_protein) + '\n')
        chisq_fh.write(','.join('%s' % x for x in TAD_nb_sig_pro) + '\n')
        chisq_fh.write(TAD_CELL + ' near boundaries of rRNA genes\n')
        chisq_fh.write(','.join('%s' % x for x in TAD_nb_rRNA) + '\n')
        chisq_fh.write(','.join('%s' % x for x in TAD_nb_sig_rRNA) + '\n')

#####################################
# Plot line graphs for each gene type
#####################################
font_size = 12
plt.rcParams['font.size'] = font_size
for g in gene_types:
    nameplot = GENE_LINE + g + '_' + TAD_CELL + '.png'
    g_loc = geneLocations[g]
    fig, ax = plt.subplots(1)
    g_loc.plot()
    plt.xlabel('Bins (Normalized TAD length)')
    plt.ylabel('Frequency')
    plt.title('Gene Location across ' + TAD_CELL + ' TADs\n' + g)
    plt.grid(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    plt.savefig(nameplot)
    plt.close()

#####################################
# Plot Chromosome Specific Line Graphs
#####################################
# Separate Graphs
for chrom in Chromosome_Specific_TAD.keys():
    if chrom in ['Boundary', 'Y']:
        continue
    nameplot = GENE_CHROM + chrom + '_' + TAD_CELL + '.png'
    chrom_loc = chromLocations[chrom]
    fig, ax = plt.subplots(1)
    chrom_loc.plot()
    plt.xlabel('Bins (Normalized TAD length)')
    plt.ylabel('Frequency')
    plt.title('Gene location (Start Site) across ' + TAD_CELL +
              ' TADs\nChromosome ' + str(chrom))
    plt.grid(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    plt.savefig(nameplot)
    plt.close()

# Combined Graph
fig, ax = plt.subplots(1)
plt.xlabel('Bins (Normalized TAD length)')
plt.ylabel('Frequency')
plt.title('Gene location (Start Site) across TADs\nby Chromosome')
plt.grid(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')

for chrom in Chromosome_Specific_TAD.keys():
    if chrom in ['Boundary']:
        continue
    chrom_loc = chromLocations[chrom]
    chrom_loc.plot()

plt.savefig(GENE_ALL_CHROM)
plt.close()

####################################
# Obtain the number of genes in each TAD
####################################
howmanygenesinTAD = []
TADlength = []
for key, value in TADdictgenes.iteritems():
    if key in ['Y', 'Boundary']:
        continue
    for tadkey, tadval in TADdictgenes[key].iteritems():
        howmanygenesinTAD.append(len(tadval))
        tadname = parse_TAD_name(tadkey)
        TADlength.append((tadname[2] - tadname[1])/1000)

####################################
# Output Histograms
####################################
# Plot number of genes in tads
font_size = 9
plt.rcParams['font.size'] = font_size
n, bins, patches = plt.hist(howmanygenesinTAD, NUM_BINS, normed=1,
                            facecolor='green', alpha=0.75)
plt.xlabel('Number of genes')
plt.ylabel('Frequency')
plt.title('Distribution of number of genes in ' + TAD_CELL + ' TADs')
plt.axis([min(howmanygenesinTAD), max(howmanygenesinTAD), 0, 0.06])
plt.grid(False)
addtext = 'mean = ' + str(round(np.mean(howmanygenesinTAD), 2)) + '\n' + \
          'std = ' + str(round(np.std(howmanygenesinTAD), 2))
plt.text(130, .055, addtext, verticalalignment='top',
         horizontalalignment='right')
plt.savefig(GENE_HIST)
plt.close()

# Plot how long TADs are
n, bins, patches = plt.hist(TADlength, NUM_BINS * 2, normed=1,
                            facecolor='green', alpha=0.75)
plt.xlabel('TAD length (kb)')
plt.ylabel('Frequency')
plt.title('How long are ' + TAD_CELL + ' TADs?')
addtext = 'mean = ' + str(round(np.mean(TADlength), 2)) + '\n' + 'std = ' + \
          str(round(np.std(TADlength), 2))
if TAD_CELL == 'hESC':
    plt.axis([min(TADlength), max(TADlength), 0, 0.0023])
    plt.text(4000, .0021, addtext, verticalalignment='top',
             horizontalalignment='right')
else:
    plt.axis([min(TADlength), max(TADlength), 0, 0.0015])
    plt.text(4000, .0011, addtext, verticalalignment='top',
             horizontalalignment='right')
plt.grid(False)
plt.savefig(TAD_LEN_HIST)
plt.close()

####################################
# PART 3 - Repeats #################
####################################

####################################
# Process repeat picklefile
####################################
all_repeat_types = []
for key in TADdictrepeats.keys():
    if key in ['Y', 'Boundary']:
        continue
    for tadkey in TADdictrepeats[key].keys():
        repeat_elements = TADdictrepeats[key][tadkey]
        repeat_elements = repeat_elements['repeat'].tolist()
        all_repeat_types = all_repeat_types + repeat_elements

all_repeat_types = set(all_repeat_types)
repeat_types = []
for repeat in all_repeat_types:
    if '?' not in repeat:
        repeat_types.append(repeat)

repeat_types.append('all')

####################################
# Observe location of repeats in TADs
####################################
# Initialize a dictionary that will store location information for repeats
TADrepeatpos = dict.fromkeys(repeat_types)
for repeatkey in repeat_types:
    # Create a list of zeros of len(NUM_BINS)
    TADrepeatpos[repeatkey] = [0] * (NUM_BINS)

# Initialize a chromosome specific dictionary holding bin specific info
Chromosome_Specific_TAD = dict.fromkeys(TADdictrepeats.keys())
for chromkey in Chromosome_Specific_TAD.keys():
    if chromkey not in ['Boundary']:
        # Create a list of zeros of len(NUM_BINS)
        Chromosome_Specific_TAD[chromkey] = [0] * (NUM_BINS)

# Obtain the repeat position information
for key, value in TADdictrepeats.iteritems():
    if key in ['Y', 'Boundary']:
        continue
    for tadkey, tadval in TADdictrepeats[key].iteritems():
        tadkey = parse_TAD_name(tadkey)
        for repeat in range(len(tadval)):
            # repeat_info stores [gene_type, chromosome, start_pos, divergence]
            repeat_info = parse_repeat_info(tadval.ix[repeat, :])
            repeat_name = repeat_info[0]

            if '?' not in repeat_name:
                binID = ID_TAD_bins(tadkey, NUM_BINS, repeat_info[2],
                                    IDtype='gene')
                if binID == -1:
                    continue

                # if the start is exactly on a TAD boundary, put it in the
                # previous bin
                if binID == 50:
                    binID = 49

                # Add to bin of specific repeat
                TADrepeatpos[repeat_name][binID] += 1

                # Add the bin number of all repeats
                TADrepeatpos['all'][binID] += 1

                # Add a count for chromsome specificity
                Chromosome_Specific_TAD[key][binID] += 1

# Convert to pandas dataframe
repeatLocations = pd.DataFrame(TADrepeatpos)
chromLocations = pd.DataFrame(Chromosome_Specific_TAD)

####################################
# Obtain the number of Repeats in each TAD
####################################
howmanyrepeatsinTAD = []
tad_names = []
for key, value in TADdictrepeats.iteritems():
    if key not in ['Y', 'Boundary']:
        for tadkey, tadval in TADdictrepeats[key].iteritems():
            howmanyrepeatsinTAD.append(len(tadval))
            tad_names.append(tadkey)

howmanyrepeatsinTAD = np.array(howmanyrepeatsinTAD)

####################################
# Output Histograms
####################################
# Plot number of genes in tads
n, bins, patches = plt.hist(howmanyrepeatsinTAD, NUM_BINS, normed=1,
                            facecolor='green', alpha=0.75)
plt.xlabel('Number of Repeat Elements')
plt.ylabel('Frequency')
plt.title('Distribution of number of Repeat Elements in ' + TAD_CELL + ' TADs')
plt.axis([min(howmanyrepeatsinTAD), max(howmanyrepeatsinTAD), 0, 0.0008])
plt.grid(False)
addtext = 'mean = ' + str(round(np.mean(howmanyrepeatsinTAD), 2)) + '\n' + \
          'std = ' + str(round(np.std(howmanyrepeatsinTAD), 2))
plt.text(5000, .00055, addtext, verticalalignment='top',
         horizontalalignment='right')
plt.savefig(REP_HIST)
plt.close()

#####################################
# Plot line graphs for each gene type
#####################################
font_size = 10
plt.rcParams['font.size'] = font_size
for r in repeat_types:
    fn = r.replace('/', '_')
    nameplot = REP_TYPE + fn + '.png'
    g_loc = repeatLocations[r]
    fig, ax = plt.subplots(1)
    g_loc.plot()
    plt.xlabel('Bins (Normalized TAD length)')
    plt.ylabel('Frequency')
    plt.title('Repeat Location across ' + TAD_CELL + ' TADs\n' + r)
    plt.grid(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    plt.savefig(nameplot)
    plt.close()

#####################################
# Plot Chromosome Specific Line Graphs
#####################################
# Separate Graphs
for chrom in Chromosome_Specific_TAD.keys():
    if chrom in ['Y', 'Boundary']:
        continue
    nameplot = REP_CHROM + chrom + '_' + TAD_CELL + '.png'
    chrom_loc = chromLocations[chrom]
    fig, ax = plt.subplots(1)
    chrom_loc.plot()
    plt.xlabel('Bins (Normalized TAD length)')
    plt.ylabel('Frequency')
    plt.title('Repeat location across ' + TAD_CELL + 'TADs\nChromosome ' +
              str(chrom))
    plt.grid(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    plt.savefig(nameplot)
    plt.close()

# Combined Graph
fig, ax = plt.subplots(1)
plt.xlabel('Bins (Normalized TAD length)')
plt.ylabel('Frequency')
plt.title('Repeat location across ' + TAD_CELL + ' TADs\nby Chromosome')
plt.grid(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')

for chrom in Chromosome_Specific_TAD.keys():
    if chrom not in ['Y', 'Boundary']:
        chrom_loc = chromLocations[chrom]
        chrom_loc.plot()

plt.savefig(REP_ALL_CHROM)
plt.close()
