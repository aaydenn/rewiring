library(dplyr)
library(JGL)
library(igraph)

load("Data/mouse-data-features.RData")
load("Data/DE_miRNAs.RData")
source("Code/tuning_parameter.r")

# modify diet column of features
features$diet <-
  ifelse(features$diet %in% c("ICRR","ICRRF"),  "ICR", features$diet)

# features$diet <-
#   ifelse(features$diet == "ICRR",  "ICR", features$diet)
# features$diet <-
#   ifelse(features$diet == "ICRRF", "ICR", features$diet)

mouse.Data |> # filter miRNAs according to DE miRNAs not in Brain
  as.data.frame() |>
  filter(rownames(mouse.Data) %in% c(DE.blood, DE.mfp)) |> t() -> arrays

# colnames(arrays) <- str_remove(colnames(arrays), "mmu-")
colnames(features) <- c("diet", "age", "tissue", "array.design")


# function for network comparison of diets for given tissue

tissues <- c(Blood = "Blood", MFP = "MFP")

diets <- c(ICR = "ICR", CCR = "CCR", AL = "AL")

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

names(res$theta) <- diets




rmat.blood <- 
  theta2rmat(cfgl$Blood$result$theta)

rmat.mfp <- 
  theta2rmat(cfgl$MFP$result$theta)

# 1="ICR", 2="CCR", 3="AL"

theta <- list(
  al.blood = rmat.blood[[3]],
  ccr.blood = rmat.blood[[2]],
  icr.blood = rmat.blood[[1]],
  al.mfp = rmat.mfp[[3]],
  ccr.mfp = rmat.mfp[[2]],
  icr.mfp = rmat.mfp[[1]]
)

theta.graphs <-
  lapply(theta, function(x) {
    G = (abs(x) > 1e-5) * 1
    diag(G) = 0
    degree = apply(G, 1, sum)
    g = G[degree > 0, degree > 0]
    graph = graph_from_adjacency_matrix(g, mode = "undirected", weighted = TRUE, diag = FALSE)
  })

# 1="ICR", 2="CCR", 3="AL"

delta <- 
  list(
    d.ccr.al.blood = rmat.blood[[2]] - rmat.blood[[3]],
    d.icr.al.blood = rmat.blood[[1]] - rmat.blood[[3]],
    d.icr.ccr.blood = rmat.blood[[1]] - rmat.blood[[2]],
    d.ccr.al.mfp = rmat.mfp[[2]] - rmat.mfp[[3]],
    d.icr.al.mfp = rmat.mfp[[1]] - rmat.mfp[[3]],
    d.icr.ccr.mfp = rmat.mfp[[1]] - rmat.mfp[[2]]
  )

delta.graphs <-
  lapply(delta, function(x) {
    G = (abs(x) > 1e-5) * 1
    diag(G) = 0
    degree = apply(G, 1, sum)
    g = G[degree > 0, degree > 0]
    graph = graph_from_adjacency_matrix(g, mode = "undirected", weighted = TRUE, diag = FALSE)
  })

centrality.table <- function(net) { # function for centrality measures for given network
  data.frame(degree = strength(net) |> signif(digits=3), 
             # closeness = closeness(net) |> signif(digits=3),
             # betweenness = betweenness(net) |> signif(digits=3),
             eigen = eigen_centrality(net)$vector |> signif(digits=3),
             knn = knn(net, V(net))$knn |> signif(digits=3)) |> 
    rownames_to_column(var = "miRNA") |> 
    mutate(hub = ifelse(knn<degree, 1,0)) |> 
    arrange(desc(eigen))
}

theta.centrality <- lapply(theta.graphs, centrality.table)
delta.centrality <- lapply(delta.graphs, centrality.table)

theta.hub <- lapply(theta.centrality, function(x) x |> filter(hub > 0) |> arrange(desc(eigen)))
delta.hub <- lapply(delta.centrality, function(x) x |> filter(hub > 0) |> arrange(desc(eigen)))

save(cfgl, theta.graphs, delta.graphs, theta.centrality, delta.centrality, 
     file = "Data/cfgl.RData")
