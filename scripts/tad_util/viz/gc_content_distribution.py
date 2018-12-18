"""
(C) 2016 Gregory Way
gc_content_distribution.py

Description:
Observe GC Content distributions and locations across TADs (hg19)

Usage:
Is called by 'ANALYSIS.sh' but can also be run through command line.

        python gc_content_distribution.py -g <GENOME>

Where <GENOME> can be either 'hg' or 'mm' for human and mouse.

Output:
GC Content distribution and histogram across TADs as a .png file
"""

import sys
sys.path.append('bin/')
import random
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import argparse
from util import parse_TAD_name
from util import ID_TAD_bins
from util import parse_repeat_info
plt.figure.max_open_warning = 0

random.seed(123)
##################
# Load Command Arguments
##################
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--genome', help='Can be either hg19 or mm9')
args = parser.parse_args()

##################
# Define Constants
##################
NUM_BINS = 50
GENOME = args.genome

if GENOME == 'hg':
    FASTA_LOC = 'data/hg/hg19_fasta/'

    TAD_DICT_A = 'index/REPEATSindex_hg19_hESC.p'
    TAD_DICT_B = 'index/REPEATSindex_hg19_IMR90.p'

    OUTPUT_FH_A = 'figures/hg19/hg19_hESC_gc_distribution.png'
    OUTPUT_FH_B = 'figures/hg19/hg19_IMR90_gc_distribution.png'

    PLT_TITLE_A = 'hg19 GC Content Across hESC TADs'
    PLT_TITLE_B = 'hg19 GC Content Across IMR90 TADs'

    OUTPUT_FH_DIV_A = 'figures/hg19/repeats/divergence/hg19_hESC_divergence_'
    OUTPUT_FH_DIV_B = 'figures/hg19/repeats/divergence/hg19_IMR90_divergence_'

    PLT_TITLE_DIV_A = 'hg19 hESC Divergence: '
    PLT_TITLE_DIV_B = 'hg19 IMR90 Divergence: '

elif GENOME == 'mm':
    FASTA_LOC = 'data/mm/mm9_fasta/'

    TAD_DICT_A = 'index/REPEATSindex_mm9_mESC.p'
    TAD_DICT_B = 'index/REPEATSindex_mm9_cortex.p'

    OUTPUT_FH_A = 'figures/mm9/mm9_mESC_gc_distribution.png'
    OUTPUT_FH_B = 'figures/mm9/mm9_cortex_gc_distribution.png'

    PLT_TITLE_A = 'mm9 GC Content Across mESC TADs'
    PLT_TITLE_B = 'mm9 GC Content Across cortex TADs'

    OUTPUT_FH_DIV_A = 'figures/mm9/repeats/divergence/mm9_mESC_divergence_'
    OUTPUT_FH_DIV_B = 'figures/mm9/repeats/divergence/mm9_cortex_divergence_'

    PLT_TITLE_DIV_A = 'mm9 mESC Divergence: '
    PLT_TITLE_DIV_B = 'mm9 Cortex Divergence: '

else:
    raise ValueError('Please input either "hg" or "mm" for -g')

####################################
# Set matplotlib defaults
####################################
fig_size = 4, 3
margin = 2.5
font_size = 10
plt.rcParams['figure.figsize'] = fig_size
plt.rcParams['font.size'] = font_size
x0, x1, y0, y1 = plt.axis()
plt.axis((x0 + margin, x1 + margin, y0 + margin, y1 + margin))

####################################
# Define Functions
####################################


def load_fasta(chrom, fasta_loc):
    """
    Retrieve fasta file

    Arguments:
    :param chrom: the chromosome of interest (format 'chr#')
    :param fasta_loc: the location that stores the hg19 fasta files

    Output:
    fasta file for the given chromosome
    """

    chrom_fa = fasta_loc + chrom + '.fa'
    record = SeqIO.read(open(chrom_fa), 'fasta')

    nucleotides = str(record.seq)
    # FASTA file has upper and lowercase letters
    # lower case = repetative elements (ignore for now)
    nucleotides = nucleotides.upper()

    return nucleotides


def parse_fasta_tad(tadname, nuc):
    """
    Retrieve nucleotide sequence of tad

    Arguments:
    :param tadname: of the format tadid:tadstart-tadend
    :param nuc: the nucleotides for the same chromosome as the TAD

    Output:
    a list of nucleotides in the TAD
    """

    tadname = parse_TAD_name(tadname)
    start = tadname[1]
    end = tadname[2]

    tadsequence = nuc[start:end]

    return tadsequence


def determine_gc_content(seq, start, end):
    """
    Determine the gc content for a given sequence given bin coordinates

    Arguments:
    :param seq: a nucleotide sequence
    :param start: where to subset the sequence
    :param end: where to subset the sequence

    Output:
    A count of GC content within the specific coordinates
    """

    if len(seq) < end:
        end = len(seq)
    subset_seq = seq[start:end]

    A = subset_seq.count('A')
    C = subset_seq.count('C')
    G = subset_seq.count('G')
    T = subset_seq.count('T')
    N = subset_seq.count('N')

    known_length = A + C + G + T

    if known_length + N == N:
        GC = 0.5
    else:
        GC = (G + C) / float(known_length)

    return float(GC)


def split_TAD_bins(tadlength, num_bins):
    """
    Return a list of coordinates to partition the TAD

    Arguments:
    :param tadlength: how long the TAD is
    :param num_bins: how many bins to split

    Output:
    a list of tuples with the locations for starting and ending bins
    """

    avgbin = tadlength / num_bins
    remainder = tadlength % num_bins
    if remainder > 0:
        randadd = random.sample(range(0, num_bins), remainder)
    else:
        randadd = []

    return_list = []
    current_idx = 0
    for binID in range(0, num_bins):
        next_idx = current_idx + avgbin
        if binID in randadd:
            return_list.append([current_idx + 1, next_idx + 1])
            current_idx = next_idx + 1
        else:
            return_list.append([current_idx + 1, next_idx])
            current_idx = next_idx
    return return_list


def get_gc_content(tads, seq, chrom, num_bins):
    """
    Determine the gc content of all TADs across TAD bins

    Arguments:
    :param tads: The lookup dictionary for TAD locations
    :param seq: the seqIO object for the chromosome fasta
    :param chrom: the chromosome we're currently using
    :param num_bins: how to split the TAD bins for determining gc content

    Output:
    A list of gc contant of all TADs across bins
    """
    for TAD in tads[chrom].keys():
        tad_sequence = parse_fasta_tad(TAD, seq)
        tad_length = len(tad_sequence)

        # Get the TAD bins
        tad_bins = split_TAD_bins(tad_length, num_bins)

        # Now, loop over the TAD bins and extract GC Content
        tad_gc = []

        for coord_idx, coord in enumerate(tad_bins):
            start = tad_bins[coord_idx][0]
            end = tad_bins[coord_idx][1]
            gc = determine_gc_content(tad_sequence, start, end)
            tad_gc.append(gc)

        return tad_gc


def identify_TAD_repeats(tads):
    """
    Return the divergence
    """
    all_repeat_types = []
    for key in tads.keys():
        if key in ['Y', 'Boundary']:
            continue
        for tadkey in tads[key].keys():
            repeat_elements = tads[key][tadkey]
            repeat_elements = repeat_elements['repeat'].tolist()
            all_repeat_types = all_repeat_types + repeat_elements

    all_repeat_types = set(all_repeat_types)
    repeat_types = []
    for repeat in all_repeat_types:
        if '?' not in repeat:
            repeat_types.append(repeat)

    repeat_types.append('all')
    # Initialize dictionary sthat will store divergence information for repeats
    TADrepeat_div = dict.fromkeys(repeat_types)
    for repeatkey in repeat_types:
        # Create a list of zeros of len(NUM_BINS)
        TADrepeat_div[repeatkey] = [[] for _ in range(NUM_BINS)]

    return TADrepeat_div


def plot_gc_content(gc_distrib, title, out_fh):
    """
    Output gc content distribution plots

    Arguments:
    :param gc_distrib: a list of gc content across TADs
    :param title: title for plot
    :param out_fh: the file handle to save

    Output:
    None -- saves plot to file
    """
    gc_distrib = pd.DataFrame(gc_distrib)

    fig, ax = plt.subplots()
    gc_distrib.boxplot(column=range(0, 50), notch=True, positions=range(0, 50),
                       return_type='dict')
    ax.set_xticklabels('' * NUM_BINS)
    ax.tick_params(axis=u'both', which=u'both', length=0)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.grid(False)
    plt.xlabel('Bins (Normalized TAD length)')
    plt.ylabel('GC Percentage')
    plt.title(title)
    plt.savefig(out_fh)
    plt.close()


def plot_divergence(div, title_base, out_base):
    """
    Output divergence plots

    Arguments:
    :param div: a list of divergence across TADs
    :param title_base: base name for title of plot
    :param out_fh: the base file handle to save

    Output:
    None -- saves plot to file
    """
    for divkey in div.keys():

        title = title_base + divkey
        out_fh = out_base + divkey.replace('/', '_') + '.png'

        fig, ax = plt.subplots()
        plt.boxplot(div[divkey])
        ax.set_xticklabels('' * NUM_BINS)
        ax.tick_params(axis=u'both', which=u'both', length=0)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.grid(False)
        plt.xlabel('Bins (Normalized TAD length)')
        plt.ylabel('Percent Divergence from Consensus Sequence')
        plt.title(title)
        plt.savefig(out_fh)
        plt.close()

####################################
# Load Data
####################################
TADrepeats_A = pd.read_pickle(TAD_DICT_A)
TADrepeats_div_A = identify_TAD_repeats(TADrepeats_A)
TADrepeats_B = pd.read_pickle(TAD_DICT_B)
TADrepeats_div_B = identify_TAD_repeats(TADrepeats_B)

####################################
# Analysis
####################################
gc_distrib_A = []
gc_distrib_B = []

for chrom in TADrepeats_A.keys():
    if chrom in ['Boundary', 'Y']:
        continue

    # Build key to lookup fasta
    lookup_chrom = 'chr' + str(chrom)

    # Get the TAD sequence:
    chr_seq = load_fasta(lookup_chrom, FASTA_LOC)

    # Determine GC content and Divergence
    for TAD, tadval in TADrepeats_A[chrom].iteritems():

        tad_sequence = parse_fasta_tad(TAD, chr_seq)
        tad_length = len(tad_sequence)

        # Get the TAD bins
        tad_bins = split_TAD_bins(tad_length, NUM_BINS)

        # Now, loop over the TAD bins and extract GC Content
        tad_gc = []

        for coord_idx, coord in enumerate(tad_bins):
            start = tad_bins[coord_idx][0]
            end = tad_bins[coord_idx][1]
            gc = determine_gc_content(tad_sequence, start, end)
            tad_gc.append(gc)

        gc_distrib_A.append(tad_gc)

        # Now, store divergence
        tadkey = parse_TAD_name(TAD)
        for repeat in range(len(tadval)):
            # repeat_info stores [gene_type, chromosome, start_pos, divergence]
            repeat_info = parse_repeat_info(tadval.ix[repeat, :])

            repeat_name = repeat_info[0]
            repeat_div = float(repeat_info[3])

            if '?' not in repeat_name:
                binID = ID_TAD_bins(tadkey, NUM_BINS, repeat_info[2],
                                    IDtype='gene')

                if binID == -1:
                    continue

                # if the start is exactly on a TAD boundary, put it in the
                # previous bin
                if binID == 50:
                    binID = 49

                # Add divergence info to specific bin
                TADrepeats_div_A[repeat_name][binID].append(repeat_div)

                # Add the bin number of all repeats
                TADrepeats_div_A['all'][binID].append(repeat_div)

    # Determine GC content and Divergence
    for TAD, tadval in TADrepeats_B[chrom].iteritems():
        tad_sequence = parse_fasta_tad(TAD, chr_seq)
        tad_length = len(tad_sequence)

        # Get the TAD bins
        tad_bins = split_TAD_bins(tad_length, NUM_BINS)

        # Loop over the TAD bins and extract GC Content
        tad_gc = []

        for coord_idx, coord in enumerate(tad_bins):
            start = tad_bins[coord_idx][0]
            end = tad_bins[coord_idx][1]
            gc = determine_gc_content(tad_sequence, start, end)
            tad_gc.append(gc)

        gc_distrib_B.append(tad_gc)

        # Store divergence
        tadkey = parse_TAD_name(TAD)
        for repeat in range(len(tadval)):
            # repeat_info stores [gene_type, chromosome, start_pos, divergence]
            repeat_info = parse_repeat_info(tadval.ix[repeat, :])
            repeat_name = repeat_info[0]
            repeat_div = repeat_info[3]
            if '?' not in repeat_name:
                binID = ID_TAD_bins(tadkey, NUM_BINS, repeat_info[2],
                                    IDtype='gene')
                if binID == -1:
                    continue

                # if the start is exactly on a TAD boundary, put it in the
                # previous bin
                if binID == 50:
                    binID = 49

                # Add divergence info to specific bin
                TADrepeats_div_B[repeat_name][binID].append(repeat_div)

                # Add the bin number of all repeats
                TADrepeats_div_B['all'][binID].append(repeat_div)

####################################
# Plot
####################################
plot_gc_content(gc_distrib_A, PLT_TITLE_A, OUTPUT_FH_A)
plot_gc_content(gc_distrib_B, PLT_TITLE_B, OUTPUT_FH_B)

plot_divergence(TADrepeats_div_A, PLT_TITLE_DIV_A, OUTPUT_FH_DIV_A)
plot_divergence(TADrepeats_div_B, PLT_TITLE_DIV_B, OUTPUT_FH_DIV_B)
