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
args <- commandArgs(trailingOnly = T)

# Parse command argumentsx
evidence_file <- args[1]
trait <- unlist(strsplit(basename(evidence_file), "_"))[1]
venn_output_file <- file.path("results", paste0("venn_", trait, ".tiff"))

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
