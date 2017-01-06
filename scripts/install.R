# 2016 Gregory Way
# TAD Pathways
# scripts/install.R
#
# The script installs all R dependencies required for the pipeline
#
# Usage:
# Run once before implementing pipeline: "Rscript scripts/install.R"

library("methods")

mirror <- "http://cran.us.r-project.org"
install.packages("checkpoint", repos = mirror)

library("checkpoint")
checkpoint("2016-02-25")

library("readr")
library("VennDiagram")
library("dplyr")
library("gridExtra")
library("optparse")

source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt", suppressUpdates = TRUE)

library("biomaRt")

sink("sessionInfo.txt")
sessionInfo()
sink()
