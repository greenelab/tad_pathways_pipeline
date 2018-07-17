#!/bin/bash
# TAD Pathways - Analytical code for TAD based analysis and visualization

# Usage:
# run_pipeline.sh usage: run_pipeline.sh [-h] [-t TAD_BOUNDARY] [-f TAD_FILE]
# optional arguments:
# -h, --help            show help message and exit
# -t, --TAD-Boundary    boundary cell type
# -f, --TAD-File        path to 3-column-tab-separated TAD domain bed file
#-s, --SNP_File         path to the snp file

# process arguments
function usage()
{ 
    echo "usage: run_pipeline.sh [-h] [-t TAD_BOUNDARY] [-f TAD_FILE]"
    echo ""
    echo "optional arguments:"
    echo "-h, --help            show this help message and exit"
    echo "-t, --TAD-Boundary    boundary cell type"
    echo "-f, --TAD-File        path to 3-column-tab-separated TAD domain bed file"
    echo "-s, --SNP_File        path to the snp file"	
}

while [ "$1" != "" ]; do
    PARAM=$1
    VALUE=$2
    case $PARAM in
        -h | --help)
            usage
            exit
            ;;
        -t | --TAD-Boundary)
            tad_cell=$VALUE
            ;;
        -f | --TAD-File)
            tad_file=$VALUE
            ;;
        -s | --SNP_File)
            candidate_snp_file=$VALUE
            ;;
        *)
            echo "ERROR: unknown parameter \"$PARAM\""
            usage
            exit 1
            ;;
    esac
    shift
    shift
done

## Define filenames
candidate_snp_location_file='results/'$tad_cell'_location.tsv'
candidate_snp_tad_file='results/'$tad_cell'_tad_results.tsv'
nearest_gene_file='results/'$tad_cell'_tad_results_nearest_gene.tsv'
trait=$tad_cell
evidence_file='results/'$tad_cell'_gene_evidence.csv'
pathway_p_values_file='gestalt/'$tad_cell'_pvals.tsv'

# Generate index files (maps to TAD identifiers to enable fast lookup)
# 1000G SNP / genes / repeat elements
python scripts/generate_index_files.py --TAD-Boundary $tad_cell  --TAD-File $tad_file

# Visualize SNPs and Genes in TADs
# Output histograms and line graphs of SNP/Gene/Repeat locations in TADs
# and gc content distribution across human and mouse tads
python scripts/visualize_genomic_elements.py --TAD-Boundary $tad_cell
python scripts/visualize_gc_and_divergence.py --TAD-Boundary $tad_cell --TAD-File $tad_file

# Map SNPs to genomic location
Rscript --vanilla scripts/build_snp_list.R \
        --snp_file $candidate_snp_file \
        --output_file $candidate_snp_location_file

# Build a customized genelist based on SNP locations
python scripts/build_custom_tad_genelist.py \
        --snp_data_file $candidate_snp_location_file \
        --output_file $candidate_snp_tad_file \
        --TAD-Boundary $tad_cell

# Perform WebGestalt pathway analysis and parse results
Rscript --vanilla scripts/webgestalt_run.R \
        --tad_genelist_file $candidate_snp_tad_file \
        --output_name $trait

# Construct an evidence file - Nearest gene to gwas or not
python scripts/construct_evidence.py \
            --trait $trait \
            --gwas $nearest_gene_file \
            --pathway $pathway_p_values_file

# Summarize the evidence file
python scripts/summarize_evidence.py \
            --evidence $evidence_file \
            --snps $candidate_snp_tad_file \
            --output_file 'results/'$tad_cell'_gene_evidence_summary.tsv'

# Visualize overlap in TAD pathways curation
R --no-save --args $evidence_file \
                < scripts/integrative_summary.R 

