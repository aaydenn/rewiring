library(dplyr)
library(JGL)
# library(igraph)

RESULT_FOLDER <- "Data"
MOUSE_DATA <- paste(RESULT_FOLDER, "mouse-data-features.RData", sep = "/")
DE_DATA    <- paste(RESULT_FOLDER, "DE_miRNAs.RData", sep = "/")

load(MOUSE_DATA)
load(DE_DATA)
source("Code/tuning_parameter.R")

colnames(features) <- c("diet", "age", "tissue", "array.design")

mouse.Data |> # filter miRNAs according to DE miRNAs not in Brain
  as.data.frame() |>
  filter(rownames(mouse.Data) %in% c(DE.blood, DE.mfp)) |> 
  t() -> arrays

tissues <- c(Blood = "Blood", MFP = "MFP")
diets   <- c(ICR = "ICR", CCR = "CCR", AL = "AL")

filter.arrays <- function(diet, tissue, features, arrays) {
  arrays[features$tissue == tissue &
           features$diet %in% c(diet, "BASELINE"), ]
}

make.array.list <- function(tissue, diets, features,  arrays) {
  return(lapply(diets, filter.arrays, tissue, features, arrays))
}

M <- lapply(tissues, make.array.list, diets, features, arrays)
allM <- c(M$Blood, M$MFP)
M <- c(Blood = M$Blood, MFP = M$MFP, both = allM)

# res <- JGL::JGL(M$Blood, lambda1 = 0.05, lambda2 = 0.25, return.whole.theta = T)
# "getJGLTuningParamResult" is from tuning_parameter.R script
ans <- lapply(M, getJGLTuningParamResult,
                                l1vec = seq(from=0.2, to=0.4, by=0.01),
                                l2vec = seq(from=0.01, to=0.03, by=0.002))

