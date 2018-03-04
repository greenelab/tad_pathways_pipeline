# 2016 Gregory Way
# TAD Pathways
# scripts/install.R
#
# The script installs all R dependencies required for the pipeline
#
# Usage:
# Run once before implementing pipeline: "Rscript scripts/install.R"

mirror <- "http://cran.us.r-project.org"

install.packages("VennDiagram", repos = mirror)
install.packages("WebGestaltR", repos = mirror)

library("readr")
library("VennDiagram")
library("dplyr")
library("tidyr")
library("gridExtra")
library("optparse")
library("WebGestaltR")
library("biomaRt")

sink("sessionInfo.txt")
sessionInfo()
sink()
