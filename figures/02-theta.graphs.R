library(igraph)
library(stringr)
library(graphlayouts)

load("result/fgl.RData")

make.palette <- function(net,centrality,col1="white",col2="darkblue"){
  palette <- colorRampPalette(c(col1,col2))
  max <- max(centrality, na.rm=TRUE)
  vector <- centrality / max
  res <- length(V(net))
  palette(res)[as.numeric(cut(vector, breaks = res))]
}

colors <- lapply(theta.graphs, function(x) make.palette(x,eigen_centrality(x)$vector))



V(theta.graphs$Blood.AL)$name <- str_remove(V(theta.graphs$Blood.AL)$name, "mmu-")
V(theta.graphs$Blood.CCR)$name <- str_remove(V(theta.graphs$Blood.CCR)$name, "mmu-")
V(theta.graphs$Blood.ICR)$name <- str_remove(V(theta.graphs$Blood.ICR)$name, "mmu-")

V(theta.graphs$Blood.AL)$color <- colors$Blood.AL
V(theta.graphs$Blood.CCR)$color <- colors$Blood.CCR
V(theta.graphs$Blood.ICR)$color <- colors$Blood.ICR

V(theta.graphs$Blood.AL)$size <- eigen_centrality(theta.graphs$Blood.AL)$vector
V(theta.graphs$Blood.CCR)$size <- eigen_centrality(theta.graphs$Blood.CCR)$vector
V(theta.graphs$Blood.ICR)$size <- eigen_centrality(theta.graphs$Blood.ICR)$vector



V(theta.graphs$MFP.AL)$name <- str_remove(V(theta.graphs$MFP.AL)$name, "mmu-")
V(theta.graphs$MFP.CCR)$name <- str_remove(V(theta.graphs$MFP.CCR)$name, "mmu-")
V(theta.graphs$MFP.ICR)$name <- str_remove(V(theta.graphs$MFP.ICR)$name, "mmu-")

V(theta.graphs$MFP.AL)$color <- colors$MFP.AL
V(theta.graphs$MFP.CCR)$color <- colors$MFP.CCR
V(theta.graphs$MFP.ICR)$color <- colors$MFP.ICR

V(theta.graphs$MFP.AL)$size <- eigen_centrality(theta.graphs$MFP.AL)$vector
V(theta.graphs$MFP.CCR)$size <- eigen_centrality(theta.graphs$MFP.CCR)$vector
V(theta.graphs$MFP.ICR)$size <- eigen_centrality(theta.graphs$MFP.ICR)$vector

# layouts <- lapply(theta.graphs, layout_with_fr)
# lay <- merge_coords(theta.graphs, layouts)
# g <- disjoint_union(theta.graphs)
# plot(g, layout = lay, 
#      vertex.size = V(g)$size*5,
#      vertex.color=V(g)$color,vertex.label.color="black")

#' ----------------------------------------
#' Blood coexpression network visualization

png(filename = "figures/theta.blood.al.png", width = 8, height = 8, units = "in", res = 320)
plot(
  theta.graphs$Blood.AL, 
  vertex.label.cex = 0.75,
  # vertex.label = NA,
  vertex.size = eigen_centrality(theta.graphs$Blood.AL)$vector*5,
  vertex.color = alpha(colors$Blood.AL,0.85), 
  vertex.label.color = "black",
  vertex.label.dist = 0.5,
  edge.width = 0.75,
  # edge.color = alpha("grey32",0.75),
  layout = as.matrix(graphlayouts::layout_igraph_stress(theta.graphs$Blood.AL)[,c(1,2)]),
  main = "AL"
)
dev.off()

png(filename = "figures/theta.blood.ccr.png", width = 8, height = 8, units = "in", res = 320)
plot(
  theta.graphs$Blood.CCR, 
  vertex.label.cex = 0.75,
  # vertex.label = NA,
  vertex.size = eigen_centrality(theta.graphs$Blood.CCR)$vector*5,
  vertex.color = alpha(colors$Blood.CCR,0.85), 
  vertex.label.color = "black",
  vertex.label.dist = 0.5,
  edge.width = 1,
  # edge.color = "grey32",
  layout = as.matrix(graphlayouts::layout_igraph_stress(theta.graphs$Blood.CCR)[,c(1,2)]),
  main = "CCR"
)
dev.off()

png(filename = "figures/theta.blood.icr.png", width = 8, height = 8, units = "in", res = 320)
plot(
  theta.graphs$Blood.ICR, 
  vertex.label.cex = 0.75,
  # vertex.label = NA,
  vertex.size = eigen_centrality(theta.graphs$Blood.ICR)$vector*5,
  vertex.color = alpha(colors$Blood.ICR,0.85), 
  vertex.label.color = "black",
  vertex.label.dist = 0.5,
  edge.width = 1,
  # edge.color = "grey32",
  layout = as.matrix(graphlayouts::layout_igraph_stress(theta.graphs$Blood.ICR)[,c(1,2)]),
  main = "ICR"
)
dev.off()


#' --------------------------------------
#' MFP coexpression network visualization

png(filename = "figures/theta.mfp.al.png", width = 8, height = 8, units = "in", res = 320)
plot(
  theta.graphs$MFP.AL, 
  vertex.label.cex = 0.75,
  # vertex.label = NA,
  vertex.size = eigen_centrality(theta.graphs$MFP.AL)$vector*5,
  vertex.color = alpha(colors$MFP.AL,0.85), 
  vertex.label.color = "black",
  vertex.label.dist = 0.5,
  edge.width = 0.75,
  # edge.color = alpha("grey32",0.75),
  layout = as.matrix(graphlayouts::layout_igraph_stress(theta.graphs$MFP.AL)[,c(1,2)]),
  main = "AL"
)
dev.off()

png(filename = "figures/theta.mfp.ccr.png", width = 8, height = 8, units = "in", res = 320)
plot(
  theta.graphs$MFP.CCR, 
  vertex.label.cex = 0.75,
  # vertex.label = NA,
  vertex.size = eigen_centrality(theta.graphs$MFP.CCR)$vector*5,
  vertex.color = alpha(colors$MFP.CCR,0.85), 
  vertex.label.color = "black",
  vertex.label.dist = 0.5,
  edge.width = 1,
  # edge.color = "grey32",
  layout = as.matrix(graphlayouts::layout_igraph_stress(theta.graphs$MFP.CCR)[,c(1,2)]),
  main = "CCR"
)
dev.off()

png(filename = "figures/theta.mfp.icr.png", width = 8, height = 8, units = "in", res = 320)
plot(
  theta.graphs$MFP.ICR, 
  vertex.label.cex = 0.75,
  # vertex.label = NA,
  vertex.size = eigen_centrality(theta.graphs$MFP.ICR)$vector*5,
  vertex.color = alpha(colors$MFP.ICR,0.85), 
  vertex.label.color = "black",
  vertex.label.dist = 0.5,
  edge.width = 1,
  # edge.color = "grey32",
  layout = as.matrix(graphlayouts::layout_igraph_stress(theta.graphs$MFP.ICR)[,c(1,2)]),
  main = "ICR"
)
dev.off()
