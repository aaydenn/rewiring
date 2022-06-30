library(dplyr)
library(igraph)

load("G:/Drive'ım/CR-miRNA/Results/2022-06-27-jgl-optimal.Rdata")

theta <- list(
    Blood.AL = ans2$jgl.res$theta$`Blood AL`,
    Blood.CCR = ans2$jgl.res$theta$`Blood CCR`,
    Blood.ICR = ans2$jgl.res$theta$`Blood ICR`,
    MFP.AL = ans2$jgl.res$theta$`MFP AL`,
    MFP.CCR = ans2$jgl.res$theta$`MFP CCR`,
    MFP.ICR = ans2$jgl.res$theta$`MFP ICR`
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
#' Function for differential adjacency matrix

diffmat <- function(theta1, theta2) {
  G1 <- (abs(theta1) > 1e-5) * 1
  G2 <- (abs(theta2) > 1e-5) * 1
  G <- abs(G2-G1)
  
  colnames(G) <- rownames(G) <- rownames(theta1)
  
  return(G)
  
}

# 1="ICR", 2="CCR", 3="AL"

delta <- list(
  D.Blood.CCR.AL = diffmat(theta$Blood.CCR, theta$Blood.AL),
  D.Blood.ICR.AL = diffmat(theta$Blood.ICR, theta$Blood.AL),
  D.Blood.ICR.CCR = diffmat(theta$Blood.ICR, theta$Blood.CCR)
  ,
  D.MFP.CCR.AL = diffmat(theta$MFP.CCR, theta$MFP.AL),
  D.MFP.ICR.AL = diffmat(theta$MFP.ICR, theta$MFP.AL),
  D.MFP.ICR.CCR = diffmat(theta$MFP.ICR, theta$MFP.AL)
  ,
  D.AL.MFP.Blood = diffmat(theta$MFP.AL,theta$Blood.AL),
  D.CCR.MFP.Blood = diffmat(theta$MFP.CCR,theta$Blood.CCR),
  D.ICR.MFP.Blood = diffmat(theta$MFP.ICR,theta$Blood.ICR)
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

topology.table <- function(net) { # function for global topology measures for given network
  data.frame(edge = length(E(net)), 
         node = length(V(net)), 
         clustering.coef = transitivity(net),
         density = edge_density(net))
}


theta.centrality <- lapply(theta.graphs, centrality.table)

openxlsx::write.xlsx(theta.centrality,file = "tables/theta.centrality.xlsx")


theta.topology <- lapply(theta.graphs, topology.table)

do.call(rbind,theta.topology) |> 
  rownames_to_column("condition") |> 
  write_delim(file = "tables/theta.topology.tsv")




delta.centrality <- lapply(delta.graphs, centrality.table)

openxlsx::write.xlsx(delta.centrality,file = "tables/delta.centrality.xlsx")


delta.topology <- lapply(delta.graphs, topology.table)

do.call(rbind,delta.topology) |> 
  rownames_to_column("condition") |> 
  write_delim(file = "tables/delta.topology.tsv")




save(theta.graphs, delta.graphs, 
     theta.centrality, delta.centrality, 
     theta.topology, delta.topology,
     file = "result/fgl.RData")
