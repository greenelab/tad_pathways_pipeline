# 2017 Gregory Way
# scripts/gene_evidence_venn.R

# Description: 
# Summarize the number of candidate genes identified given 3 different
# aggregation methods - 1) Nearest genes 2) LD window genes 3) TAD genes

# Usage:
#   Rscript --vanilla scripts/gene_evidence_venn.R

# Output:
# Venn diagram of BMD gene aggregation methods summary

# File names
bmd_snps_file <- file.path("data", "gwas_tad_snps",
                           "Bone_mineral_density_hg19_SNPs.tsv")
tad_genes_file <- file.path("data", "gwas_tad_genes",
                            "Bone_mineral_density_hg19_SNPs_TAD_genelists.tsv")
ld_genes_file <- file.path("data", "BMD_LDwindow_genes.tsv")

# Load data and obtain genelists
bmd_snps_df <- readr::read_tsv(bmd_snps_file)
bmd_snps <- unique(bmd_snps_df$rs)
# length(bmd_snps) = 70

# Obtain all nearest genes:
# intergenic SNPs are mapped to both upstream and downstream nearest genes
bmd_genes <- bmd_snps_df$gene
bmd_nearest_genes <- c()
for (gene in bmd_genes) {
  gene <- unlist(strsplit(gene, ","))
  bmd_nearest_genes <- c(bmd_nearest_genes, gene)
}
bmd_nearest_genes <- unique(bmd_nearest_genes)

# TAD based genelist
tad_genes_df <- readr::read_tsv(tad_genes_file)
tad_genes_df <- tad_genes_df[tad_genes_df$gene_type == "protein_coding", ]
tad_genes_df <- tad_genes_df[tad_genes_df$db == "HAVANA", ]
tad_genes <- unique(tad_genes_df$gene_name)

# LD based genelist
ld_genes_df <- readr::read_tsv(ld_genes_file)
ld_genes_df <- ld_genes_df[ld_genes_df$gene_type == "protein_coding", ]
ld_genes_df <- ld_genes_df[ld_genes_df$db == "HAVANA", ]
ld_genes <- unique(ld_genes_df$gene_name)

# Make a Venn diagram of intersecting genes
candidate_genes <- unique(c(ld_genes, tad_genes, bmd_nearest_genes))

bmd_nearest_genes_idx <- match(bmd_nearest_genes, candidate_genes)
tad_genes_idx <- match(tad_genes, candidate_genes)
ld_genes_idx <- match(ld_genes, candidate_genes)

venn_list <- list()
venn_list[["Nearest"]] <- bmd_nearest_genes_idx
venn_list[["TAD"]] <- tad_genes_idx
venn_list[["LD"]] <- ld_genes_idx

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
output_venn <- VennDiagram::venn.diagram(venn_list, filename = NULL,
                                         col = "transparent",
                                         fill = c("Red", "Yellow", "Blue"),
                                         cat.cex = 1.5, euler.d = TRUE, cex = 2)

dir.create("figures")
gene_agg_filename <- file.path("figures", "gene_aggregation_venn_diagram.pdf")
pdf(gene_agg_filename, height = 6, width = 6)
grid::grid.draw(output_venn)
dev.off()

# Pathway overlap
get_path_genes <- function(file_name, pathway = "skeletal system development") {
  # Given a filename and webgestalt output, compile pathway genes
  path_df <- readr::read_tsv(file_name)
  path_df <- path_df[path_df$go_name == pathway, ]
  path_genes <- unique(path_df$symbol)
  return(path_genes)
}

tad_path_file <- file.path("gestalt", "bmd_complete_gestalt.tsv")
ld_path_file <- file.path("gestalt", "bmd_LD_complete_gestalt.tsv")
near_path_file <- file.path("gestalt", "bmd_nearest_gene_complete_gestalt.tsv")

tad_path_genes <- get_path_genes(tad_path_file)
ld_path_genes <- get_path_genes(ld_path_file)
near_path_genes <- get_path_genes(near_path_file)

# Make a Venn diagram of intersecting genes
path_genes <- unique(c(ld_path_genes, near_path_genes, tad_path_genes))

bmd_nearest_genes_idx <- match(near_path_genes, path_genes)
bmd_ld_genes_idx <- match(ld_path_genes, path_genes)
bmd_tad_genes_idx <- match(tad_path_genes, path_genes)

venn_list <- list()
venn_list[["Near Pathway"]] <- bmd_nearest_genes_idx
venn_list[["TAD_Pathway"]] <- bmd_tad_genes_idx
venn_list[["LD Pathway"]] <- bmd_ld_genes_idx

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
output_venn <- VennDiagram::venn.diagram(venn_list, filename = NULL,
                                         col = "transparent",
                                         fill = c("Red", "Yellow", "Blue"),
                                         cat.cex = 1.5, euler.d = TRUE, cex = 2)
path_filename <- file.path("figures", "pathway_aggregation_venn_diagram.pdf")
pdf(path_filename, height = 6, width = 6)
grid::grid.draw(output_venn)
dev.off()
