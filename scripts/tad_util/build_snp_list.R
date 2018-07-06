# 2016 Gregory Way
# scripts/tad_util/build_snp_list.R

# Description: 
# Maps an input file of SNPs (rsIDs) to their respective genomic coordinates
# based on GRCh37 and dbSNP141 and the ENSEMBL biomaRt bioconductor package.
# The input SNP file may have multiple columns but each column MUST have a
# column name. The input SNP file also must be comma separated (.csv).

# Usage:
# The script is run in the command line

#           Rscript --vanilla scripts/tad_util/build_snp_list.R

# and takes two positional arguments as input:

#         --snp_file        <LOCATION_OF_SNP.csv>
#         --output_file     <MAPPING_FILE.tsv>

# Output:
# A tab separated file of rsIDs mapped to genomic coordinates

library(checkpoint)
suppressMessages(checkpoint("2016-02-25"))

library(dplyr)

# Load in command arguments
option_list <- list(optparse::make_option(c("-s", "--snp_file"),
                                          type = "character",
                                          help = "Location of rsID SNP file"),
                    optparse::make_option(c("-o", "--output_file"),
                                          type = "character",
                                          help = "File name to save results"))

opt_parser <- optparse::OptionParser(option_list = option_list);
opt <- optparse::parse_args(opt_parser);

# Load arguments
snp_file <- opt$snp_file
output_file <- opt$output_file

# Map rsids to genomic location
snps <- readr::read_csv(snp_file)
snp <- biomaRt::useEnsembl(biomart = "snp", GRCh = 37,
                           dataset = "hsapiens_snp")

results_full <- list()
for (group in 1:ncol(snps)) {
  results <- biomaRt::getBM(attributes = c("refsnp_id", "refsnp_source",
                                           "chr_name", "chrom_start",
                                           "chrom_end", "minor_allele"),
                            filters = "snp_filter", values = snps[, group],
                            mart = snp)
  results <- results %>% dplyr::filter(!duplicated(refsnp_id))
  results$Group <- colnames(snps)[group]
  results_full[[group]] <- results
}

# Process results and write to file
results_full <- do.call(rbind, results_full)
results_full$chr_name <- as.numeric(results_full$chr_name)
results_full <- results_full %>% dplyr::filter(!is.na(chr_name))
colnames(results_full) <- c("snp", "refsnp_source", "chrom", "position",
                            "snp_end", "minor_allele", "group")
write.table(results_full, output_file, row.names = FALSE, sep = "\t")
