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

# Using conda version 4.4.11
conda activate tad_pathways
```

Now, a `TAD_Pathways` analysis can proceed. Follow an example pipeline to work
from an existing GWAS or the custom pipeline example for insight on how to run
`TAD_Pathways` on user curated SNPs.
 
### Examples

We provide an example for a TAD pathways analysis pipeline. To run this example:

```bash
source activate tad_pathways

# Example using custom input SNPs
bash example_pipeline_custom.sh
```

### General Usage

To perform a `TAD_Pathways` analysis, uses need to spicify 3 inputs:

1. name of the tad cell:

   E.g.: 'hESC'

2. path to the TAD domain file:

   The TAD domain file is a 3-column tab-separated bed file. The first column is the chromsome number. The second column is the  start position of the tad. And the third position is the end position of the tad.

   E.g.: [`hESC_domains_hg19.bed`](hESC_domains_hg19.bed)

3. path to the SNPs file

   The SNPs file is a comma separated text file. The first row of the text file should have group names and
subsequent rows should list the rs numbers of interest. There can be manycolumns with variable length rows.

   E.g.: [`custom_example.csv`](custom_example.csv)

| Group 1 | Group 2 |
| ------- | ------- |
| rs12345 | rs67891 |
| rs19876 | rs54321 |
| ...     | ...     |

Then, perform the following steps:

```bash
source activate tad_pathways

bash run_pipeline.sh --TAD-Boundary hESC \
                     --TAD-File hESC_domains_hg19.bed \
                     --SNP-File custom_example.csv
```

The output of these steps are Group specific text files with all genes in TADs
harboring an input SNP. See
[`example_pipeline_custom.sh`](example_pipeline_custom.sh) for more details.

### Contact

For all questions and bug reporting please file a
[GitHub issue](https://github.com/greenelab/tad_pathways/issues)

For all other questions contact Casey Greene at csgreene@mail.med.upenn.edu or
Struan Grant at grants@email.chop.edu
