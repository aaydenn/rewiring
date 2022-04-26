#!/usr/bin/env Rscript

# libFolder <- paste(getwd(), "/lib", sep = "")
# outdir <- dir.create(libFolder)

if (!require("devtools", character.only = TRUE)) {
      install.packages("devtools", dependencies = TRUE)
      library("devtools", character.only = TRUE)
}

cran.pkgs <- c(
    "tidyverse", "DiffGraph", 
    "igraph", "pbapply", 
    "optparse", "randomcoloR"
)

bioc.pkgs <- c(
    "oligo", "limma"
)

install_cran(cran.pkgs)
install_bioc(bioc.pkgs)
install_github("soumyabrataghosh/pd.mirna.4.1")


