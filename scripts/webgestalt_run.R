# 2017 Gregory Way
# scripts/webgestalt_run.R
#
# Description:
# Runs WebGestalt Pathway Analysis on Input genelist
#
# Usage:
# The script is run in the command line
#
#           Rscript --vanilla scripts/webgestalt_run.R
#
# and takes five arguments as input:
#
#         --tad_genelist_file       <LOCATION_OF_GWAS_FILE.csv>
#         --output_name             <FILENAME_ABBREVIATION>
#         --pathway                 <STRING OF PATHWAY DATABASE>
#         --topten                  <IF INCLUDED, OUTPUT TOP 10 ENRICHED SETS>
#         --output_directory        <DIRECTORY TO SAVE RESULTS>
#
# Output:
# Significantly overrepresented pathways from a WebGestalt Analysis

library(WebGestaltR)
library(tidyr)
library(dplyr)

# Load in command arguments
option_list <- list(
  optparse::make_option(c("-t", "--tad_genelist_file"),
                        type = "character",
                        help = "Location of TAD genes"),
  optparse::make_option(c("-o", "--output_name"),
                        type = "character",
                        help = "Abbrev. to save results"),
  optparse::make_option(c("-p", "--pathway"),
                        type = "character",
                        default = "geneontology_Biological_Process",
                        help = "pathway to consider"),
  optparse::make_option(c("-u", "--topten"),
                        action = "store_true",
                        help = "output top pathways without significance cut"),
  optparse::make_option(c("-d", "--output_directory"),
                        type = "character",
                        help = "Directory to save results",
                        default = "gestalt")
)

opt_parser <- optparse::OptionParser(option_list = option_list);
opt <- optparse::parse_args(opt_parser);

# Load arguments
tad_gene_file <- opt$tad_genelist_file
output_name <- opt$output_name
pathway <- opt$pathway
topten <- opt$topten
output_dir <- opt$output_directory

output_pval_file <- file.path(output_dir, paste0(output_name, "_pvals.tsv"))
output_path_file <- file.path(output_dir, paste0(output_name, "_gestalt.tsv"))

gene_df <- readr::read_tsv(tad_gene_file)
genes <- gene_df$gene_name

if (topten) {
  sigMethod <- "top"
} else {
  sigMethod <- "fdr"
}

webgestalt_output <- WebGestaltR(enrichMethod = "ORA",
                                 enrichDatabase = pathway,
                                 organism = "hsapiens",
                                 interestGene = genes,
                                 interestGeneType = "genesymbol",
                                 minNum = 4,
                                 sigMethod = sigMethod,
                                 fdrMethod = "BH",
                                 is.output = TRUE,
                                 outputDirectory = output_dir,
                                 referenceSet = "genome",
                                 projectName = output_name)

# Process output files
if (!is.null(webgestalt_output)) {
  webgestalt_output <- webgestalt_output %>%
    tidyr::separate_rows(overlapGene, OverlapGene_UserID, sep = ";")
  p_val <- webgestalt_output %>% dplyr::select(description, PValue, FDR)
  p_val <- p_val[!duplicated(p_val), ]

  colnames(webgestalt_output) <- c("go_id", "go_name", "link", "count",
                                   "observed", "expected", "R", "pval",
                                   "adjP", "overlapGene", "symbol")
  colnames(p_val) <- c("go_name", "pval", "adjP")

  write.table(p_val, output_pval_file, sep = "\t", row.names = FALSE)
  write.table(webgestalt_output, output_path_file, sep = "\t", row.names = FALSE)

}
