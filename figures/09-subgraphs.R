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


mir374 <- "mmu-miR-374c-3p"


subgraph(theta.graphs$Blood.AL, unlist(neighborhood(theta.graphs$Blood.AL, nodes = mir374))) |> 
  myplot("cadetblue4")

subgraph(theta.graphs$Blood.CCR, unlist(neighborhood(theta.graphs$Blood.CCR, nodes = mir374))) |> 
  myplot("darkolivegreen4")

subgraph(theta.graphs$Blood.ICR, unlist(neighborhood(theta.graphs$Blood.ICR, nodes = mir374))) |> 
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




fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}


svg("figures/subgraphs_manuscript.svg",8,8)

par(mfrow=c(3,3))


# blood miR-505-5p

subgraph(theta.graphs$Blood.AL, unlist(neighborhood(theta.graphs$Blood.AL, nodes = "mmu-miR-505-5p"))) |> myplot("cadetblue4")

subgraph(theta.graphs$Blood.CCR, unlist(neighborhood(theta.graphs$Blood.CCR, nodes = "mmu-miR-505-5p"))) |> myplot("darkolivegreen4")

subgraph(theta.graphs$Blood.ICR, unlist(neighborhood(theta.graphs$Blood.ICR, nodes = "mmu-miR-505-5p"))) |> myplot("coral4")



# mfp miR-23a-3p

subgraph(theta.graphs$MFP.AL, unlist(neighborhood(theta.graphs$MFP.AL, nodes = "mmu-miR-23a-3p"))) |> myplot("cadetblue4")

subgraph(theta.graphs$MFP.CCR, unlist(neighborhood(theta.graphs$MFP.CCR, nodes = "mmu-miR-23a-3p"))) |> myplot("darkolivegreen4")

subgraph(theta.graphs$MFP.ICR, unlist(neighborhood(theta.graphs$MFP.ICR, nodes = "mmu-miR-23a-3p"))) |> myplot("coral4")



# mfp miR-23a-3p

subgraph(theta.graphs$Blood.AL, unlist(neighborhood(theta.graphs$Blood.AL, nodes = "mmu-miR-144-3p"))) |> myplot("cadetblue4")

subgraph(theta.graphs$Blood.CCR, unlist(neighborhood(theta.graphs$Blood.CCR, nodes = "mmu-miR-144-3p"))) |> myplot("darkolivegreen4")

graph_from_literal("mmu-miR-144-3p") |> myplot("coral4")


dev.off()
