#!/usr/bin/env Rscript
# This code should be executed in the base folder

# TODO: simplify, since we no longer use Brain data


# Input: 
#  + intermittent_sdrf.tsv: description of the experiments
#  + *.CEL : microarray results

DATA_FOLDER <- "Data"
MASTER_TABLE <- paste(DATA_FOLDER, "intermittent_sdrf_new.tsv", sep = "/")
ANNOTATION_FILE <- paste(DATA_FOLDER,
                          "miRNA-4_1-st-v1.annotations.20160922.csv.zip")

# Output:
#  + norm-data-features.RData: containing
#    + norm.Data
#    + features
#    + mouse_4.1

RESULT_FOLDER <- "Data"
FULL_DATA  <- paste(RESULT_FOLDER, "full.Data.RData", sep = "/")
NORM_DATA  <- paste(RESULT_FOLDER, "norm-data-features.RData", sep = "/")
MOUSE_DATA <- paste(RESULT_FOLDER, "mouse-data-features.RData", sep = "/")

library(oligo, warn.conflicts = FALSE)
library(pd.mirna.4.0, warn.conflicts = FALSE)
library(pd.mirna.4.1, warn.conflicts = FALSE)
library(readr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(stringr, warn.conflicts = FALSE)

# we read the experiment description
## read master table that contains all information
## Columns are …
ff <- read.table(MASTER_TABLE, sep = "\t", header = TRUE)

## Keep only cases that did not develop breast cancer
ff |> filter(disease..breast.cancer. == FALSE) -> ff
ff$array.data.path <- paste(DATA_FOLDER, ff$organism.part,
               ff$array.data.file, sep = "/")

# PROBLEM: some of the samples were measured using a different platform
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
                 exprs = cbind(assayData(raw.Data)$exprs,
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
annot_4.1 <- read_csv(ANNOTATION_FILE, comment = "#", show_col_types = FALSE)

ff |> filter(array.design=="miRNA_4_1") |>
  select(diet, age..week., organism.part, array.design) -> features1
ff |> filter(array.design=="miRNA_4_0") |>
  select(diet, age..week., organism.part, array.design) -> features0

features <- rbind(features1, features0)
colnames(features) <- c("diet", "age", "tissue")

# we merge ICRR and ICRRF 
# we ignore one week refeed of mice
features$diet[features$diet=="ICRR"]  <- "ICR"
features$diet[features$diet=="ICRRF"] <- "ICR"
features$diet <- as.factor(features$diet)
features$diet <- relevel(features$diet, ref="BASELINE")

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

# this is a table to translate probe_id to miRNA id
# TODO: data is already in `mouse_4.1`. This can be simplified
## library(miRBaseConverter)
## miRNAid <- miRNA_AccessionToName(mouse_4.1$Accession,
##                                  targetVersion = "v22")
## no_name <- is.na(miRNAid$TargetName)
## miRNAid$TargetName[no_name] <- miRNAid$Accession[no_name]
## rownames(mouse.Data) <- str_remove(miRNAid$TargetName, "mmu-") 
rownames(mouse.Data) <- mouse_4.1$`Transcript ID(Array Design)`


# save data for further use
save(mouse.Data, features, mouse_4.1, file = MOUSE_DATA)
save(norm.Data, features, annot_4.1, file = NORM_DATA)
save(full.Data, features, eset, file = FULL_DATA)
