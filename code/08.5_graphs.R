library(igraph)
library(randomcoloR)
load("Data/d.fgl.RData")

# set.seed(42)
find_communities <- function(diet, tissue) {
  d.net <- d.fgl[[diet]][[tissue]]$d.fgl$Delta.graph.connected
  d.com <- cluster_leiden(d.net, resolution_parameter = 0.05)
  return(d.com)
}

communities <-
  lapply(c(CCR = "CCR", ICR = "ICR"), function(diet)
    lapply(c(Blood = "Blood", MFP = "MFP"),
           function(tissue)
             find_communities(diet, tissue)))

draw_net_diet_tissue <- function(diet, tissue) {
  d.net <- d.fgl[[diet]][[tissue]]$d.fgl$Delta.graph.connected
  d.com <- communities[[diet]][[tissue]]
  communityColors <-
    distinctColorPalette(max(d.com$membership))[membership(d.com)]

  png(paste(diet, "AL" , tissue, "net.png", sep = "."))
  plot(d.com, d.net, mark.groups = NULL, vertex.label.cex = 0.5,
    vertex.size = degree(d.net) * 1.5,
    col = communityColors, layout = layout_with_kk,
    main = paste("DCE network of", diet, "in", tissue)
  )
  dev.off()
}

for(diet in c("CCR","ICR")) {
  for(tissue in c("Blood", "MFP")) {
    draw_net_diet_tissue(diet, tissue)
  }
}

save(communities, file = "Data/networks-communities.Rdata")

#####
#####

V(d.fgl$CCR$Blood$d.fgl$Delta.graph.connected)$name <- 
  str_remove(V(d.fgl$CCR$Blood$d.fgl$Delta.graph.connected)$name, "mmu-")
V(d.fgl$CCR$Blood$d.fgl$Delta.graph.connected)$eigen <- 
  eigen_centrality(d.fgl$CCR$Blood$d.fgl$Delta.graph.connected)$vector

res <- length(V(d.fgl$CCR$Blood$d.fgl$Delta.graph.connected))
palette <- colorRampPalette(c('white','darkblue'))
max <- max(V(d.fgl$CCR$Blood$d.fgl$Delta.graph.connected)$eigen, na.rm=TRUE)
vector <- V(d.fgl$CCR$Blood$d.fgl$Delta.graph.connected)$eigen / max
colors <- palette(res)[as.numeric(cut(vector, breaks=res))]

write_graph(graph = d.fgl$CCR$Blood$d.fgl$Delta.graph.connected,
            file = "CCR.Blood.graphml",
            format = "graphml")

png("rewiring_ccr_blood.png", 7, 7, units = "in", res = 320)
plot(
  d.fgl$CCR$Blood$d.fgl$Delta.graph.connected, 
  vertex.label.cex = 0.5,
  vertex.size = degree(d.fgl$CCR$Blood$d.fgl$Delta.graph.connected) * 1.5,
  vertex.color = alpha(colors,0.75), 
  vertex.label.color = "black",
  vertex.label.dist = 0.5,
  edge.width = 2.25,
  edge.color = "grey42",
  layout = layout_with_kk,
  main = "Rewiring network of CCR in Blood"
)
dev.off()


V(d.fgl$ICR$Blood$d.fgl$Delta.graph.connected)$name <- 
  str_remove(V(d.fgl$ICR$Blood$d.fgl$Delta.graph.connected)$name, "mmu-")
V(d.fgl$ICR$Blood$d.fgl$Delta.graph.connected)$eigen <- 
  eigen_centrality(d.fgl$ICR$Blood$d.fgl$Delta.graph.connected)$vector

res <- length(V(d.fgl$ICR$Blood$d.fgl$Delta.graph.connected))
palette <- colorRampPalette(c('white','darkblue'))
max <- max(V(d.fgl$ICR$Blood$d.fgl$Delta.graph.connected)$eigen, na.rm=TRUE)
vector <- V(d.fgl$ICR$Blood$d.fgl$Delta.graph.connected)$eigen / max
colors <- palette(res)[as.numeric(cut(vector, breaks=res))]

write_graph(graph = d.fgl$ICR$Blood$d.fgl$Delta.graph.connected,
            file = "ICR.Blood.graphml",
            format = "graphml")

png("rewiring_icr_blood.png", 7, 7, units = "in", res = 320)
plot(
  d.fgl$ICR$Blood$d.fgl$Delta.graph.connected, 
  vertex.label.cex = 0.5,
  vertex.size = degree(d.fgl$ICR$Blood$d.fgl$Delta.graph.connected) * 1.5,
  vertex.color = alpha(colors,0.75), 
  vertex.label.color = "black",
  vertex.label.dist = 0.5,
  edge.width = 2.25,
  edge.color = "grey42",
  layout = layout_with_kk,
  main = "Rewiring network of ICR in Blood"
)
dev.off()


V(d.fgl$CCR$MFP$d.fgl$Delta.graph.connected)$name <- 
  str_remove(V(d.fgl$CCR$MFP$d.fgl$Delta.graph.connected)$name, "mmu-")
V(d.fgl$CCR$MFP$d.fgl$Delta.graph.connected)$eigen <- 
  eigen_centrality(d.fgl$CCR$MFP$d.fgl$Delta.graph.connected)$vector

res <- length(V(d.fgl$CCR$MFP$d.fgl$Delta.graph.connected))
palette <- colorRampPalette(c('white','darkblue'))
max <- max(V(d.fgl$CCR$MFP$d.fgl$Delta.graph.connected)$eigen, na.rm=TRUE)
vector <- V(d.fgl$CCR$MFP$d.fgl$Delta.graph.connected)$eigen / max
colors <- palette(res)[as.numeric(cut(vector, breaks=res))]

write_graph(graph = d.fgl$CCR$MFP$d.fgl$Delta.graph.connected,
            file = "CCR.MFP.graphml",
            format = "graphml")

png("rewiring_ccr_mfp.png", 7, 7, units = "in", res = 320)
plot(
  d.fgl$CCR$MFP$d.fgl$Delta.graph.connected, 
  vertex.label.cex = 0.5,
  vertex.size = degree(d.fgl$CCR$MFP$d.fgl$Delta.graph.connected) * 1.5,
  vertex.color = alpha(colors,0.75), 
  vertex.label.color = "black",
  vertex.label.dist = 0.5,
  edge.width = 2.25,
  edge.color = "grey42",
  layout = layout_with_kk,
  main = "Rewiring network of CCR in MFP"
)
dev.off()



V(d.fgl$ICR$MFP$d.fgl$Delta.graph.connected)$name <- 
  str_remove(V(d.fgl$ICR$MFP$d.fgl$Delta.graph.connected)$name, "mmu-")
V(d.fgl$ICR$MFP$d.fgl$Delta.graph.connected)$eigen <- 
  eigen_centrality(d.fgl$ICR$MFP$d.fgl$Delta.graph.connected)$vector

res <- length(V(d.fgl$ICR$MFP$d.fgl$Delta.graph.connected))
palette <- colorRampPalette(c('white','darkblue'))
max <- max(V(d.fgl$ICR$MFP$d.fgl$Delta.graph.connected)$eigen, na.rm=TRUE)
vector <- V(d.fgl$ICR$MFP$d.fgl$Delta.graph.connected)$eigen / max
colors <- palette(res)[as.numeric(cut(vector, breaks=res))]

write_graph(graph = d.fgl$ICR$MFP$d.fgl$Delta.graph.connected,
            file = "ICR.MFP.graphml",
            format = "graphml")

png("rewiring_icr_mfp.png", 7, 7, units = "in", res = 320)
plot(
  d.fgl$ICR$MFP$d.fgl$Delta.graph.connected, 
  vertex.label.cex = 0.5,
  vertex.size = degree(d.fgl$ICR$MFP$d.fgl$Delta.graph.connected) * 1.5,
  vertex.color = alpha(colors,0.75), 
  vertex.label.color = "black",
  vertex.label.dist = 0.5,
  edge.width = 2.25,
  edge.color = "grey42",
  layout = layout_with_kk,
  main = "Rewiring network of ICR in MFP"
)
dev.off()
