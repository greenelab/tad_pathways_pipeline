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
checkpoint("2017-05-22")

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
