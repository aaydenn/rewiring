library(dplyr)
library(igraph)

theta <- 
  list(
    Blood.AL = ans$jgl.res$theta$Blood$AL,
    Blood.CCR = ans$jgl.res$theta$Blood$CCR,
    Blood.ICR = ans$jgl.res$theta$Blood$ICR,
    MFP.AL = ans$jgl.res$theta$MFP$AL,
    MFP.CCR = ans$jgl.res$theta$MFP$CCR,
    MFP.ICR = ans$jgl.res$theta$MFP$ICR
  )

make.graph <- function(x) {
  G <- (abs(x) > 1e-5) * 1
  diag(G) <- 0
  degree <- apply(G, 1, sum)
  g <- G[degree > 0, degree > 0]
  graph <-
    graph_from_adjacency_matrix(g,
                                mode = "undirected",
                                weighted = TRUE,
                                diag = FALSE)
  return(graph)
}

theta.graphs <- lapply(theta, make.graph)

#' ------------------------------------------
#' Function for differential adjecency matrix

diffnet <- function(theta1, theta2) {
  
  m <- n <- ncol(theta1)
  
  mat <- matrix(0, m, n)
  
  for (i in 1:m) {
    for (j in 1:n) {
      if (theta1[i, j] == theta2[i, j]) {
        mat[i, j] = 0
      }
      if (abs(theta1[i, j]) > 0 &
          abs(theta2[i, j]) > 0) {
        mat[i, j] = 0
      }
      if (abs(theta2[i, j]) > 0 &
          theta1[i, j] == 0) {
        mat[i, j] = 1
      }
      if (theta2[i, j] == 0 &
          abs(theta1[i, j]) > 0) {
        mat[i, j] = 1
      }
    }
  }
  
  colnames(mat) <- rownames(mat) <- rownames(theta1)
  
  return(mat)
  
}

# 1="ICR", 2="CCR", 3="AL"

delta <- 
  list(
    D.Blood.CCR.AL = diffnet(theta$Blood.CCR, theta$Blood.AL),
    D.Blood.ICR.AL = diffnet(theta$Blood.ICR, theta$Blood.AL),
    D.Blood.ICR.CCR = diffnet(theta$Blood.ICR, theta$Blood.CCR)
  )


delta.graphs <- lapply(delta, make.graph)

centrality.table <- function(net) { # function for centrality measures for given network
  data.frame(degree = strength(net) |> signif(digits=3), 
             # closeness = closeness(net) |> signif(digits=3),
             # betweenness = betweenness(net) |> signif(digits=3),
             eigen = eigen_centrality(net)$vector |> signif(digits=3),
             knn = knn(net, V(net))$knn |> signif(digits=3)) |> 
    rownames_to_column(var = "miRNA") |> 
    mutate(hub = ifelse(knn < degree, 1,0)) |> 
    arrange(desc(eigen))
}

theta.centrality <- lapply(theta.graphs, centrality.table)
delta.centrality <- lapply(delta.graphs, centrality.table)

theta.hub <- lapply(theta.centrality, function(x) x |> filter(hub > 0) |> arrange(desc(eigen)))
delta.hub <- lapply(delta.centrality, function(x) x |> filter(hub > 0) |> arrange(desc(eigen)))

save(theta.graphs, delta.graphs, 
     theta.centrality, delta.centrality, 
     theta.hub, delta.hub,
     file = "Data/cfgl.RData")
