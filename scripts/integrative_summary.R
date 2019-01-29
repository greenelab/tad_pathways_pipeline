# 2016 Gregory Way
# scripts/integrative_summary.R

# Description:
# Take as input the genes identified by the TAD pathway analysis and the
# nearest gene GWAS to determine evidence overlaps

# Usage:
# The script is run by scripts/run_pipeline with evidence file command arg

#   R --no-save --args <EVIDENCE_FILE> < scripts/integrative_summary.R

# Output:
# Venn diagram of trait specific evidence overlaps

# Load in command arguments
option_list <- list(
  optparse::make_option(c("-e", "--evidence_file"),
                        type = "character",
                        help = "TAD gene evidence"),
  optparse::make_option(c("-t", "--trait"),
                        type = "character",
                        help = "the trait that we're focusing on"),
  optparse::make_option(c("-d", "--output_directory"),
                        type = "character",
                        help = "Directory to save results",
                        default = "results")
)

opt_parser <- optparse::OptionParser(option_list = option_list);
opt <- optparse::parse_args(opt_parser);

# Parse command arguments
evidence_file <- opt$evidence_file
trait <- opt$trait
output_dir <- opt$output_directory

# Create directory
dir.create(output_dir, showWarnings = FALSE)

# Create output file name
venn_output_file <- file.path(output_dir, paste0("venn_", trait, ".tiff"))

# Read in Data
evidence_genes <- readr::read_csv(evidence_file)

# Prepare for VennDiagram input
tad_genes <- c()
eqtl_genes <- c()
gwas_genes <- c()
for (gene in 1:nrow(evidence_genes)) {
  assignment <- evidence_genes$evidence[gene]
  if (assignment == "tad") {
    tad_genes <- c(tad_genes, gene)
  } else if (assignment == "gwas") {
    gwas_genes <- c(gwas_genes, gene)
  } else if (assignment == "gwas_tad") {
    tad_genes <- c(tad_genes, gene)
    gwas_genes <- c(gwas_genes, gene)
  }
}

venn_list <- list("GWAS" = gwas_genes, "TAD Pathway" = tad_genes)

# Output Venn Diagram
VennDiagram::venn.diagram(x = venn_list,
                          filename = venn_output_file,
                          fill = c("red", "blue"),
                          height = 1500,
                          width = 1500,
                          scaled = FALSE,
                          cat.pos = 10)
