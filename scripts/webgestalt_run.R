# 2017 Gregory Way
# scripts/webgestalt_run.R

# Description: 
# Runs WebGestalt Pathway Analysis on Input genelist

# Usage:
# The script is run in the command line

#           Rscript --vanilla scripts/webgestalt_run.R

# and takes two positional arguments as input:

#         --tad_genelist_file       <LOCATION_OF_GWAS_FILE.csv>
#         --output_name             <FILENAME_ABBREVIATION>

# Output:
# A tab separated file of rsIDs mapped to genomic coordinates

library(checkpoint)
suppressMessages(checkpoint("2017-05-22"))

library(WebGestaltR)
library(tidyr)
library(dplyr)

# Load in command arguments
option_list <- list(optparse::make_option(c("-t", "--tad_genelist_file"),
                                          type = "character",
                                          help = "Location of TAD genes"),
                    optparse::make_option(c("-o", "--output_name"),
                                          type = "character",
                                          help = "Abbrev. to save results"))

opt_parser <- optparse::OptionParser(option_list = option_list);
opt <- optparse::parse_args(opt_parser);

# Load arguments
tad_gene_file <- opt$tad_genelist_file
output_name <- opt$output_name

output_pval_file <- file.path("gestalt", paste0(output_name, "_pvals.tsv"))
output_path_file <- file.path("gestalt", paste0(output_name, "_gestalt.tsv"))

gene_df <- readr::read_tsv(tad_gene_file)
genes <- gene_df$gene_name

webgestalt_output <- WebGestaltR(enrichMethod = "ORA",
                                 organism = "hsapiens",
                                 interestGene = genes,
                                 interestGeneType = "genesymbol",
                                 minNum = 4,
                                 fdrMethod = "BH",
                                 is.output = TRUE,
                                 outputDirectory = "gestalt",
                                 referenceSet = "genome",
                                 projectName = output_name)

# Process output files
webgestalt_output <- webgestalt_output %>%
  tidyr::separate_rows(overlapGene, OverlapGene_UserID, sep = ";")
p_val <- webgestalt_output %>% dplyr::select(description, PValue) 
p_val <- p_val[!duplicated(p_val), ]

colnames(webgestalt_output) <- c("go_id", "go_name", "link", "count",
                                 "observed", "expected", "R", "pval",
                                 "adjP", "overlapGene", "symbol")
colnames(p_val) <- c("go_name", "adjP")

write.table(p_val, output_pval_file, sep = "\t", row.names = FALSE)
write.table(webgestalt_output, output_path_file, sep = "\t", row.names = FALSE)
