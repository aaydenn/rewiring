library(circlize)
load("data/go_kegg_list.RData")
load("result/targets.Rdata")

# theta.top.go2mirna <- lapply(theta.go.top, function(x) lapply(GO.list.BP[x$GO_Term], function(y) unlist(symbol2mir_id[names(symbol2mir_id) %in% y])))

mychord <- function(tab) {
  
  circos.par(gap.degree = 2, cell.padding = c(0, 0, 0, 0))
  
  chordDiagram(
    tab,
    annotationTrack = "grid",
    transparency = 0.75,
    # link.lwd = 2,
    # link.border = "gray",
    # link.zindex = rank(tab),
    preAllocateTracks = list(track.height = 0.2)
  )
  
  circos.track(
    track.index = 1,
    panel.fun = function(x, y) {
      circos.text(
        CELL_META$xcenter,
        CELL_META$ylim[1],
        CELL_META$sector.index,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5)
      )
    },
    bg.border = NA
  )
  
  circos.clear()
}



theta.go.top <- lapply(theta.go, function(x) x |> head(10L))

delta.go.top <- lapply(delta.go, function(x) x |> head(10L))


tab1 <- GO.table.BP[GO.table.BP$GO_Term %in% theta.go.top$Blood.AL$GO_Term,] |> mutate(symbol = Mouse_Marker_Symbol)

mir2go <- targets |> right_join(tab1, by = "symbol") |> select(mirna,GO_Term) |> distinct() |> na.omit()
       
mychord(mir2go)


# filter targets for 98 mirna
targets <- targets.unique[targets.unique$mirna %in% colnames(theta[[1]]),]

# filter top targeted genes ( quantile(0.99) = 4.31 )
targets.top <- targets |> group_by(symbol) |> summarise(n=n()) |> arrange(desc(n)) |> mutate(z=(n-mean(n))/sd(n)) |> filter(z>4.31)


# data frame of top targeted genes and their targetşng miRNAs
df <- targets[targets$symbol %in% targets.top$symbol,c("mirna","symbol")] |> mutate(mirna=str_remove(mirna,pattern = "mmu-"))


png(filename = "figures/target.genes.png",width = 8,height = 8,units = "in",res = 300)

circos.par(gap.degree = 1, cell.padding = c(0, 0, 0, 0))
chordDiagram(
  df,
  annotationTrack = "grid",
  transparency = 0.75,
  preAllocateTracks = list(track.height = 0.2)
)
circos.track(
  track.index = 1,
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1],
      CELL_META$sector.index,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5)
    )
  },
  bg.border = NA
)
circos.clear()

dev.off()

