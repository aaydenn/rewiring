theta <- list(
  al.blood = rmat.blood[[3]],
  ccr.blood = rmat.blood[[2]],
  icr.blood = rmat.blood[[1]],
  al.mfp = rmat.mfp[[3]],
  ccr.mfp = rmat.mfp[[2]],
  icr.mfp = rmat.mfp[[1]]
)

make.graph <- function(x) {
    G <- (abs(x) > 1e-5) * 1
    diag(G) <- 0
    degree <- apply(G, 1, sum)
    g <- G[degree > 0, degree > 0]
    graph <- graph_from_adjacency_matrix(g, mode = "undirected", weighted = TRUE, diag = FALSE)
    return(graph)
  }
theta.graphs <- lapply(theta, make.graph)

# 1="ICR", 2="CCR", 3="AL"

delta <- list(
    d.ccr.al.blood = rmat.blood[[2]] - rmat.blood[[3]],
    d.icr.al.blood = rmat.blood[[1]] - rmat.blood[[3]],
    d.icr.ccr.blood = rmat.blood[[1]] - rmat.blood[[2]],
    d.ccr.al.mfp = rmat.mfp[[2]] - rmat.mfp[[3]],
    d.icr.al.mfp = rmat.mfp[[1]] - rmat.mfp[[3]],
    d.icr.ccr.mfp = rmat.mfp[[1]] - rmat.mfp[[2]]
  )

delta.graphs <- lapply(delta, function(x) {
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
