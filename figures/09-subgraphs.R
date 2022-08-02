library(igraph)

load("result/targets.Rdata")
load("result/theta_delta.RData")



# igraph object of target data
mirna_target <- graph_from_data_frame(targets.unique[c("mirna","symbol")],directed = FALSE)


# union of theta and target
target_theta <- lapply(theta.graphs, function(x) x %u% mirna_target) |> 
  
  # set type of nodes (bipartite: mirna and gene)
  lapply(function(x) set_vertex_attr(x,"type",value = ifelse(grepl("mmu-", V(x)$name), "mirna", "gene")))


# subgraph of selected mirnas
mirlist <- list(
  `mmu-miR-374c-3p`="mmu-miR-374c-3p", `mmu-miR-505-5p`="mmu-miR-505-5p", 
  `mmu-let-7e-5p`="mmu-let-7e-5p", `mmu-miR-6991-5p`="mmu-miR-6991-5p", 
  `mmu-miR-144-3p`="mmu-miR-144-3p",`mmu-miR-668-3p` = "mmu-miR-668-3p",
  `mmu-miR-6939-5p` = "mmu-miR-6939-5p",`mmu-miR-466h-3p` = "mmu-miR-466h-3p",
  `mmu-miR-324-3p` = "mmu-miR-324-3p",`mmu-miR-331-3p` = "mmu-miR-331-3p",
  `mmu-miR-7044-5p` = "mmu-miR-7044-5p"
)


target_theta_sub <- lapply(target_theta, function(x) 
  lapply(mirlist, function(y) 
    subgraph(x, vids = unlist(neighborhood(x, nodes = y)))
    ))

myplot <- function(graph, color) {
  plot(graph, layout = layout_with_dh, vertex.label = str_remove(V(graph)$name, "mmu-"),vertex.label.color = "black", vertex.frame.color = NA, vertex.color = color, edge.color = "gray51",  edge.width = 2, curved = TRUE, margin = -0.05)
}



#' --------------------------
#' "mmu-miR-374c-3p" in blood



png("figures/mmu-miR-374c-3p.png", 8,4,"in",res = 300)

par(mfrow = c(1,3))


subgraph(theta.graphs$Blood.AL, unlist(neighborhood(theta.graphs$Blood.AL, nodes = "mmu-miR-374c-3p"))) |> 
  myplot("cadetblue4")

subgraph(theta.graphs$Blood.CCR, unlist(neighborhood(theta.graphs$Blood.CCR, nodes = "mmu-miR-374c-3p"))) |> 
  myplot("darkolivegreen4")

subgraph(theta.graphs$Blood.ICR, unlist(neighborhood(theta.graphs$Blood.ICR, nodes = "mmu-miR-374c-3p"))) |> 
  myplot("coral4")


dev.off()





#' --------------------------
#' "mmu-miR-144-3p" in blood



png("figures/mmu-miR-144-3p.png", 8,4,"in",res = 300)

par(mfrow = c(1,3))



subgraph(theta.graphs$Blood.AL, unlist(neighborhood(theta.graphs$Blood.AL, nodes = "mmu-miR-144-3p"))) |> myplot("cadetblue4")

subgraph(theta.graphs$Blood.CCR, unlist(neighborhood(theta.graphs$Blood.CCR, nodes = "mmu-miR-144-3p"))) |> myplot("darkolivegreen4")

graph_from_literal("mmu-miR-144-3p") |> myplot("coral4")

dev.off()





#' --------------------------
#' "mmu-miR-505-5p" in blood



png("figures/mmu-miR-505-5p.png", 8,4,"in",res = 300)

par(mfrow = c(1,3))



subgraph(theta.graphs$Blood.AL, unlist(neighborhood(theta.graphs$Blood.AL, nodes = "mmu-miR-505-5p"))) |> myplot("cadetblue4")

subgraph(theta.graphs$Blood.CCR, unlist(neighborhood(theta.graphs$Blood.CCR, nodes = "mmu-miR-505-5p"))) |> myplot("darkolivegreen4")

subgraph(theta.graphs$Blood.ICR, unlist(neighborhood(theta.graphs$Blood.ICR, nodes = "mmu-miR-505-5p"))) |> myplot("coral4")

dev.off()





#' --------------------------
#' "mmu-miR-23a-3p" in mfp



png("figures/mmu-miR-23a-3p.png", 8,4,"in",res = 300)

par(mfrow = c(1,3))



subgraph(theta.graphs$MFP.AL, unlist(neighborhood(theta.graphs$MFP.AL, nodes = "mmu-miR-23a-3p"))) |> myplot("cadetblue4")

subgraph(theta.graphs$MFP.CCR, unlist(neighborhood(theta.graphs$MFP.CCR, nodes = "mmu-miR-23a-3p"))) |> myplot("darkolivegreen4")

subgraph(theta.graphs$MFP.ICR, unlist(neighborhood(theta.graphs$MFP.ICR, nodes = "mmu-miR-23a-3p"))) |> myplot("coral4")

dev.off()











