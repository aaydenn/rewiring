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

layouts <- lapply(theta.graphs, layout_with_fr)

lay <- merge_coords(theta.graphs, layouts)

g <- disjoint_union(theta.graphs)

plot(g,
     layout = lay, vertex.size = V(g)$size*5,
     vertex.color=V(g)$color,vertex.label.color="black")


# color palette





png("figures/theta.blood.png", units = "in", res = 320)

par(mfrow=c(3,1), mar = c(1, 1, 1, 1))

plot(
  theta.graphs$Blood.AL, 
  vertex.label.cex = 0,
  vertex.label = NA,
  vertex.size = eigen_centrality(theta.graphs$Blood.AL)$vector*5,
  vertex.color = alpha(colors$Blood.AL,0.85), 
  vertex.label.color = "black",
  vertex.label.dist = 0.5,
  edge.width = 1,
  # edge.color = "grey32",
  layout = layout_with_fr(theta.graphs$Blood.AL),
  main = "AL"
)


plot(
  theta.graphs$Blood.CCR, 
  vertex.label.cex = 0,
  vertex.label = NA,
  vertex.size = eigen_centrality(theta.graphs$Blood.CCR)$vector*5,
  vertex.color = alpha(colors$Blood.CCR,0.85), 
  vertex.label.color = "black",
  vertex.label.dist = 0.5,
  edge.width = 1,
  # edge.color = "grey32",
  layout = layout_with_fr(theta.graphs$Blood.CCR),
  main = "CCR"
)

plot(
  theta.graphs$Blood.ICR, 
  vertex.label.cex = 0,
  vertex.label = NA,
  vertex.size = eigen_centrality(theta.graphs$Blood.ICR)$vector*5,
  vertex.color = alpha(colors$Blood.ICR,0.85), 
  vertex.label.color = "black",
  vertex.label.dist = 0.5,
  edge.width = 1,
  # edge.color = "grey32",
  layout = layout_with_fr(theta.graphs$Blood.ICR),
  main = "ICR"
)

dev.off()
