library(circlize)

load("result/targets.Rdata")
load("result/DE_miRNAs.RData")

# theta.top.go2mirna <- lapply(theta.go.top, function(x) lapply(GO.list.BP[x$GO_Term], function(y) unlist(symbol2mir_id[names(symbol2mir_id) %in% y])))




# Blood
# filter targets for 63 mirna

targets_blood <- targets.unique |> 
  filter(mirna %in% DE.blood)


# filter top targeted genes ( quantile(0.99) = 4.31 )
targets.top <- targets_blood |> 
  group_by(symbol) |> 
  summarise(n=n()) |> 
  arrange(desc(n)) |> 
  mutate(z=(n-mean(n))/sd(n)) |> 
  filter(z>quantile(z,0.95))


# data frame of top targeted genes and their targeted miRNAs
blood <- targets_blood |> filter(symbol %in% targets.top$symbol) |> 
  select(c("mirna","symbol")) |> 
  mutate(mirna=str_remove(mirna,pattern = "mmu-"))


png(filename = "figures/target.genes.blood.png", width = 8, height = 8, units = "in",res = 300)

circos.par(gap.degree = 1, cell.padding = c(0, 0, 0, 0))

chordDiagram(
  blood,
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
      adj = c(0, 0.5),
      cex = 0.75
    )
  },
  bg.border = NA
)

circos.clear()

dev.off()




# MFP
# filter targets for 35 mirna

targets_mfp <- targets.unique |> 
  filter(mirna %in% DE.mfp)



# filter top targeted genes ( quantile(0.99) = 4.31 )
targets.top <- targets_mfp |> 
  group_by(symbol) |> 
  summarise(n=n()) |> 
  arrange(desc(n)) |> 
  mutate(z=(n-mean(n))/sd(n)) |> 
  filter(z>quantile(z,0.95))


# data frame of top targeted genes and their targeted miRNAs
mfp <- targets_mfp |> filter(symbol %in% targets.top$symbol) |> 
  select(c("mirna","symbol")) |> 
  mutate(mirna=str_remove(mirna,pattern = "mmu-"))



png(filename = "figures/target.genes.mfp.png", width = 8, height = 8, units = "in",res = 300)

circos.par(gap.degree = 1, cell.padding = c(0, 0, 0, 0))

chordDiagram(
  mfp,
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
      adj = c(0, 0.05),
      cex = 0.75
    )
  },
  bg.border = NA
)

circos.clear()

dev.off()
