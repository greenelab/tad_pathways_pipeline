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
                        default = FALSE,
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

# Create folder
dir.create(output_dir)

output_pval_file <- file.path(output_dir, paste0(output_name, "_pvals.tsv"))
output_path_file <- file.path(output_dir, paste0(output_name, "_gestalt.tsv"))

gene_df <- readr::read_tsv(tad_gene_file)
genes <- gene_df$gene_name

if (topten) {
  sigMethod <- "top"
} else {
  sigMethod <- "fdr"
}

webgestalt_df <- WebGestaltR(enrichMethod = "ORA",
                             organism = "hsapiens",
                             enrichDatabase = pathway,
                             interestGene = genes,
                             interestGeneType = "genesymbol",
                             minNum = 4,
                             sigMethod = sigMethod,
                             fdrMethod = "BH",
                             isOutput = TRUE,
                             outputDirectory = output_dir,
                             referenceSet = "genome",
                             referenceGeneType = "genename",
                             projectName = output_name)

# Load ID mapping
id_map_file <- file.path(output_dir,
                         paste0("Project_", output_name),
                         paste0("interestingID_mappingTable_", output_name, ".txt"))
id_df <- readr::read_tsv(id_map_file,
                         col_types = readr::cols(.default = readr::col_character(),
                                                 entrezgene = readr::col_character()))

# Process output files
if (!is.null(webgestalt_df)) {
  webgestalt_df <- webgestalt_df %>%
    tidyr::separate_rows(overlapId, sep = ";")

  p_val <- webgestalt_df %>% dplyr::select(description, pValue, FDR)
  p_val <- p_val[!duplicated(p_val), ]

  colnames(p_val) <- c("id", "pval", "adjP")

  webgestalt_df <- webgestalt_df %>%
      dplyr::left_join(id_df, by = c("overlapId" = "entrezgene"))

  colnames(webgestalt_df) <- c("id", "term", "link", "count",
                               "observed", "expected", "R", "pval",
                               "adjP", "entrezgene", "symbol",
                               "genename", "genelink")

  readr::write_tsv(p_val, output_pval_file)
  readr::write_tsv(webgestalt_df, output_path_file)
}
