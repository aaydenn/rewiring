library(dplyr)
library(JGL)
# library(igraph)

load("Data/mouse-data-features.RData")
load("Data/DE_miRNAs.RData")
source("Code/tuning_parameter.R")


# colnames(arrays) <- str_remove(colnames(arrays), "mmu-")
colnames(features) <- c("diet", "age", "tissue", "array.design")

# modify diet column of features
features$diet <-
  ifelse(features$diet %in% c("ICRR","ICRRF"),  
         "ICR", 
         features$diet)

mouse.Data |> # filter miRNAs according to DE miRNAs not in Brain
  as.data.frame() |>
  filter(rownames(mouse.Data) %in% c(DE.blood, DE.mfp)) |> 
  t() -> arrays


# function for network comparison of diets for given tissue

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

# res <- JGL::JGL(M$Blood, lambda1 = 0.05, lambda2 = 0.25, return.whole.theta = T)
# "getJGLTuningParamResult" is from tuning_parameter.R script
ans <- lapply(M, getJGLTuningParamResult)
ans2 <- getJGLTuningParamResult(allM, l1vec = seq(from=0.2, to=0.4, by=0.01),
                                l2vec = seq(from=0.01, to=0.03, by=0.002))

