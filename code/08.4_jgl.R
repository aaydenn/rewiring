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


#' Derive the Joint Graphical Lasso estimation result by AIC tuning parameter selection
#'
#' @param GauList a list of gaussian distributed matrices
#' @param l1vec a vector of candidate values for lambda1, if NULL, the values will be 1:5/20
#' @param l2vec a vector of candidate values for lambda2, if NULL, the values will be 1:5/50
#' @return aic.vec: a table of AIC values ; jgl.res: the JGL result object
#' @export
getJGLTuningParallel <- function(GauList, l1vec = NULL, l2vec = NULL) {
  # Gaulist: samples by genes
  if (is.null(l1vec)){
    l1vec = 1:5/20
  }
  if (is.null(l2vec)){
    l2vec = 1:5/50
  }
  L1 <- rep(l1vec, each  = length(l2vec))
  L2 <- rep(l2vec, times = length(l1vec))
  jgl.res <- lapply(1:length(L1), function(i) {
    AIC_select(mat.list.t.gau = GauList, lam1 = L1[i], lam2 = L2[i],
              returnJGL = T))
  })
  return(jgl.res)
}


# res <- JGL::JGL(M$Blood, lambda1 = 0.05, lambda2 = 0.25, return.whole.theta = T)
# "getJGLTuningParamResult" is from tuning_parameter.R script
ans <- lapply(M, getJGLTuningParamResult,
                                l1vec = seq(from=0.2, to=0.4, by=0.01),
                                l2vec = seq(from=0.01, to=0.03, by=0.002))

