# TAD_Pathways

## Leveraging TADs to identify candidate genes at GWAS signals

**Gregory P. Way and Casey S. Greene - 2017**

### Summary

The repository contains data and instructions to implement a "TAD_Pathways"
analysis for over 300 different trait/disease GWAS or custom SNP lists.
TAD_Pathways uses the principles of topologically association domains (TADs) to
define where an association signal (typically a GWAS signal) can most likely
impact gene function. We use TAD boundaries as defined by
[Dixon et al. 2012](https://doi.org/10.1038/nature11082) and
[hg19 Gencode genes](ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/)
to identify which genes may be implicated. We then input this list into a
[WebGestalt Pathways Analysis](http://www.webgestalt.org/) to output
significantly associated pathways implicated by the input TAD-defined geneset.


For more specific details about the method, refer to our
[preprint](https://doi.org/10.1101/087718).

### Setup

Before you begin, download the necessary TAD based index files and
GWAS curation files and setup python environment:

```bash
bash initialize.sh

source activate tad_pathways
```

### Usage

There are two ways to implement a TAD_Pathways analysis:

1. GWAS
2. Custom

#### GWAS

Browse the `data/gwas/` directory to select a GWAS file. Each file in this
directory is a flat text file that includes a single column of genes that are
located within TADs implicated by a significant GWAS association. For complete
information on how these lists were constructed, refer to
https://github.com/greenelab/tad_pathways

Input this gene list directly into a
[WebGestalt Pathway Analysis](http://http://www.webgestalt.org/). Skip now to
the [WebGestalt step](#webgestalt-pathway-analysis).

#### Custom

Create a comma separated file where the first row of each column names the list
of snps below in subsequent rows. There can be many columns with variable
length rows.

E.g.: `custom_example.csv`

| Group 1 | Group 2 |
| ------- | ------- |
| rs12345 | rs67891 |
| rs19876 | rs54321 |
| ...     | ...     |

Then, perform the following steps:

```bash
# Map custom SNPs to genomic locations
Rscript --vanilla scripts/build_snp_list.R \
        --snp_file "custom_example.csv" \
        --output_file "mapped_results.tsv"

# Build TAD based genelists for each group
python scripts/build_custom_TAD_genelist.py \
       --snp_data_file "mapped_results.tsv" \
       --output_file "custom_tad_genelist.tsv"
```

Skip now to the the [WebGestalt step](#webgestalt-pathway-analysis).

### WebGestalt Pathway Analysis

Insert either the GWAS curated genelist or a column from the custom genelist
with the following parameters:

| Parameter | Input |
| --------- | ----- |
| Select gene ID type | *hsapiens__gene_symbol* |
| Enrichment Analysis | *GO Analysis* |
| GO Slim Classification | Yes |
| Reference Set | *hsapiens__genome* |
| Statistical Method | *Hypergeometric* |
| Multiple Test Adjustment | *BH* |
| Significance Level | *Top10* |
| Minimum Number of Genes for a Category | *4* |

Once the analysis is complete, click `Export TSV Only` and save the file as
`data/gestalt/<INSERT_TRAIT_HERE>_gestalt.tsv`. 

### Curation

Clean and tidy the output files and summarize into convenient lists of
candidate genes. These genes may or may not be the nearest gene to the GWAS
signal and all require experimental validation.

```bash
# Process WebGestalt Output
# `INSERT_TRAIT_HERE` is a string indicating the name of the file
# saved in `data/gestalt/<INSERT_TRAIT_HERE>_gestalt.tsv`
python scripts/parse_gestalt.py --trait 'INSERT_TRAIT_HERE' --process

# Output evidence tables
python scripts/construct_evidence.py \
        --trait 'BMD' \
        --eqtl 'data/eqtl/eqtl_BMD_genelist.tsv' \
        --genelist 'data/gwas_catalog/Bone_mineral_density_hg19.tsv' \
        --pathway 'skeletal system development'

# Summarize evidence
python scripts/assign_evidence_to_TADs.py \
        --evidence 'tad_pathway/BMD_gene_evidence.csv' \
        --snps 'data/gwas_TAD_location/Bone_mineral_density_hg19_SNPs.tsv' \
        --output_file 'tad_pathway/BMD_evidence_summary.tsv'

# Output venn diagrams and gene lists
R --no-save --args 'tad_pathway/BMD_gene_evidence.csv' \
        'BMD' < scripts/integrative_summary.R
```

### Examples

We implemented two examples as a proof of concept. To view these examples,
initialize the repository and run:

```bash
# Bone Mineral Density GWAS
bash bmd_example_pipeline.sh

# Type 2 Diabetes GWAS
bash t2d_example_pipeline.sh
```

