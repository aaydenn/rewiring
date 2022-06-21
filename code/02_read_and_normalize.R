#!/usr/bin/env Rscript
# This code should be executed in the Data folder
#
# Input: 
#  + intermittent_sdrf.tsv: description of the experiments
#  + *.CEL : microarray results
#
# Output:
#  + norm-data-features.RData: containing
#    + norm.Data
#    + features
#    + mouse_4.1

library(oligo, warn.conflicts = FALSE)
library(pd.mirna.4.0, warn.conflicts = FALSE)
library(pd.mirna.4.1, warn.conflicts = FALSE)
library(readr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(stringr, warn.conflicts = FALSE)

# we read the experiment description
ff <- read.table("Data/intermittent_sdrf_new.tsv", sep="\t", header=TRUE) #read master table that contains all information
ff |> filter(disease..breast.cancer. == FALSE) -> ff
ff$array.data.path <- paste("Data", ff$organism.part, ff$array.data.file, 
                            sep = "/")
# first reading the 4.1 arrays
ff |>
  filter(array.design=="miRNA_4_1") |>
  select(array.data.path) |>
  read.celfiles(pkgname = "pd.mirna.4.1") -> raw.Data
# then reading the 4.0 arrays
ff |>
  filter(array.design=="miRNA_4_0") |>
  select(array.data.path) |>
  read.celfiles(pkgname = "pd.mirna.4.0") -> raw.Data0
# merge to raw intensities into one expression
a1 <- assayData(raw.Data)$exprs
a2 <- assayData(raw.Data0)$exprs

full.Data <- new("ExpressionFeatureSet",
                 exprs=cbind(assayData(raw.Data)$exprs,
                          assayData(raw.Data0)$exprs))
  
pheno <- rbind(raw.Data@phenoData@data, raw.Data0@phenoData@data)
pheno$index <- 1:nrow(pheno)
full.Data@phenoData@data <- pheno
full.Data@annotation <- raw.Data@annotation
full.Data@protocolData@data$exprs <- c(
  raw.Data@protocolData@data$exprs,
  raw.Data0@protocolData@data$exprs)
full.Data@protocolData@data$dates <- c(
  raw.Data@protocolData@data$dates,
  raw.Data0@protocolData@data$dates)

eset <- rma(full.Data)
norm.Data <- exprs(eset)
annot_4.1 <- read_csv("Data/miRNA-4_1-st-v1.annotations.20160922.csv.zip",
                      comment = "#", show_col_types = FALSE)

ff |>
  filter(array.design=="miRNA_4_1") |>
  select(diet, age..week., organism.part, array.design) -> features1
ff |>
  filter(array.design=="miRNA_4_0") |>
  select(diet, age..week., organism.part, array.design) -> features0

features <- rbind(features1, features0)
colnames(features) <- c("diet", "age", "tissue")

annot_4.1 |> filter(`Species Scientific Name`=="Mus musculus", 
                    str_starts(Accession,"MIMAT")) -> mouse_4.1
# we keep only mouse probes
mouse.Data <- norm.Data[mouse_4.1$`Probe Set Name`,]
colnames(mouse.Data) |> 
  stringr::str_remove(".CEL") |> 
  stringr::str_remove(".cel") |> 
  stringr::str_remove(".ga") -> colnames(mouse.Data)
# change some names for consistency
colnames(mouse.Data)[colnames(mouse.Data) =="T71Blood"] <- "ICRR_T71_Blood_w80"
colnames(mouse.Data)[colnames(mouse.Data) =="T25W80Blood"] <- "ICRRF_T25_Blood_w80"
# save data for further use
save(mouse.Data, features, mouse_4.1, file = "Data/mouse-data-features.RData")
save(norm.Data, features, annot_4.1, file = "Data/norm-data-features.RData")
save(full.Data, features, eset, file = "Data/full.Data.RData")
