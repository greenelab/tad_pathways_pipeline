# TAD_Pathways

## Leveraging TADs to identify candidate genes at GWAS signals

**Gregory P. Way, Casey S. Greene, and Struan F.A. Grant - 2017**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.254190.svg)](https://doi.org/10.5281/zenodo.254190)

### Summary

The repository contains data and instructions to implement a "TAD_Pathways"
analysis for over 300 different trait/disease GWAS or custom SNP lists.

TAD_Pathways uses the principles of topologically association domains (TADs) to
define where an association signal (typically a GWAS signal) can most likely
impact gene function. We use TAD boundaries as defined by
[Dixon et al. 2012](https://doi.org/10.1038/nature11082) and
[hg19 Gencode genes](ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/)
to identify which genes may be implicated. We then perform an overrepresentation
pathway analysis to identify significantly associated pathways implicated by the
input TAD-defined geneset.

For more specific details about our method, refer to our
[short report](https://doi.org/10.1038/ejhg.2017.108 "Implicating candidate genes at GWAS signals by leveraging topologically associating domains")
at the European Journal of Human Genetics.

We also present a 6 minute video introducing the method and discussing the
experimental validation at
[EJHG-tube](http://www.nature.com/ejhg/videos/index.html).

### Setup

First, clone the repository and navigate into the top directory:

```bash
git clone git@github.com:greenelab/tad_pathways_pipeline.git
cd tad_pathways_pipeline
```

Before you begin, download the necessary TAD based index files and GWAS
curation files and setup python environment:

```bash
bash initialize.sh

conda env create --force --file environment.yml
source activate tad_pathways
```

Now, a Tad_pathways analysis can proceed. Follow an example pipeline to work
from existing GWAS or the custom pipeline example for insight on how to run
Tad_pathways on user curated SNPs.
 
### Examples

We provide three different examples for a TAD pathways analysis pipeline. To run
each of the analyses:

```bash
# Example using Bone Mineral Density GWAS
bash example_pipeline_bmd.sh

# Example using Type 2 Diabetes GWAS
bash example_pipeline_t2d.sh

# Example using custom input SNPs
bash example_pipeline_custom.sh
```

### General Usage

There are two ways to implement a TAD_Pathways analysis:

1. GWAS
2. Custom

#### GWAS

Browse the `data/gwas_tad_genes/` directory to select a GWAS file. Each file in
this directory is a tab separated text file that includes information regarding
each gene located within a signal TAD. The column `gene_name` is the
comprehensive list of all implicated genes. For complete information on how
these lists were constructed, refer to
https://github.com/greenelab/tad_pathways. 

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

The output of these steps are Group specific text files with all genes in TADs
harboring an input SNP. See
[`example_pipeline_custom.sh`](example_pipeline_custom.sh) for more details.

### Contact

For all questions and bug reporting please file a
[GitHub issue](https://github.com/greenelab/tad_pathways/issues)

For all other questions contact Casey Greene at csgreene@mail.med.upenn.edu or
Struan Grant at grants@email.chop.edu
