library(igraph)
library(graphlayouts)

load("result/diff.rank.eigen.RData")
load("result/fgl.RData")


make.palette <- function(net,centrality,col1="white",col2="goldenrod"){
  palette <- colorRampPalette(c(col1,col2))
  max <- max(centrality, na.rm=TRUE)
  vector <- centrality / max
  res <- length(V(net))
  palette(res)[as.numeric(cut(vector, breaks = res))]
}



#' ------------
#' Blood subgraph

D.Blood.CCR.AL.subgraph <- induced_subgraph(delta.graphs$D.Blood.CCR.AL, diff.rank.eigen$diff.rank.Blood.CCR.AL$miRNA)

colors <- make.palette(D.Blood.CCR.AL.subgraph, eigen_centrality(D.Blood.CCR.AL.subgraph)$vector)

png(filename = "figures/D.Blood.CCR.AL.subgraph.png", width = 8, height = 8, units = "in", res = 300)
plot(
  D.Blood.CCR.AL.subgraph,
  vertex.labe.cex = 0.75,
  vertex.label = str_remove(V(D.Blood.CCR.AL.subgraph)$name,"mmu-"),
  vertex.size = degree(D.Blood.CCR.AL.subgraph)+10,
  vertex.color = alpha(colors, 0.85),
  vertex.label.color = "black",
  layout = layout_with_fr,
  main = "CCR.AL"
)
dev.off()

D.Blood.ICR.AL.subgraph <- induced_subgraph(delta.graphs$D.Blood.ICR.AL, diff.rank.eigen$diff.rank.Blood.ICR.AL$miRNA)

colors <- make.palette(D.Blood.ICR.AL.subgraph, eigen_centrality(D.Blood.ICR.AL.subgraph)$vector)


png(filename = "figures/D.Blood.ICR.AL.subgraph.png", width = 8, height = 8, units = "in", res = 300)
plot(
  D.Blood.ICR.AL.subgraph,
  vertex.labe.cex = 0.75,
  vertex.label = str_remove(V(D.Blood.ICR.AL.subgraph)$name,"mmu-"),
  vertex.size = degree(D.Blood.ICR.AL.subgraph)+10,
  vertex.color = alpha(colors, 0.85),
  vertex.label.color = "black",
  layout = layout_with_fr,
  main = "ICR.AL"
)
dev.off()


#' ------------
#' MFP subgraph

D.MFP.CCR.AL.subgraph <- induced_subgraph(delta.graphs$D.MFP.CCR.AL, diff.rank.eigen$diff.rank.MFP.CCR.AL$miRNA)

colors <- make.palette(D.MFP.CCR.AL.subgraph, eigen_centrality(D.MFP.CCR.AL.subgraph)$vector)

png(filename = "figures/D.MFP.CCR.AL.subgraph.png", width = 8, height = 8, units = "in", res = 300)
plot(
  D.MFP.CCR.AL.subgraph,
  vertex.labe.cex = 0.75,
  vertex.label = str_remove(V(D.MFP.CCR.AL.subgraph)$name,"mmu-"),
  vertex.size = degree(D.MFP.CCR.AL.subgraph)+10,
  vertex.color = alpha(colors, 0.85),
  vertex.label.color = "black",
  layout = layout_with_fr,
  main = "CCR.AL"
)
dev.off()



D.MFP.ICR.AL.subgraph <- induced_subgraph(delta.graphs$D.MFP.ICR.AL, diff.rank.eigen$diff.rank.MFP.ICR.AL$miRNA)

colors <- make.palette(D.MFP.ICR.AL.subgraph, eigen_centrality(D.MFP.ICR.AL.subgraph)$vector)


png(filename = "figures/D.MFP.ICR.AL.subgraph.png", width = 8, height = 8, units = "in", res = 300)
plot(
  D.MFP.ICR.AL.subgraph,
  vertex.labe.cex = 0.75,
  vertex.label = str_remove(V(D.MFP.ICR.AL.subgraph)$name,"mmu-"),
  vertex.size = degree(D.MFP.ICR.AL.subgraph)+10,
  vertex.color = alpha(colors, 0.85),
  vertex.label.color = "black",
  layout = layout_with_fr,
  main = "ICR.AL"
)
dev.off()

