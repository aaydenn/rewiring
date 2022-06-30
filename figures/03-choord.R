library(circlize)
load("data/go_kegg_list.RData")

mychord <- function(tab) {
  
  circos.par(gap.degree = 1, cell.padding = c(0, 0, 0, 0))
  
  chordDiagram(
    tab,
    annotationTrack = "grid",
    transparency = 0.75,
    # link.lwd = 2,
    # link.border = "gray",
    link.zindex = rank(tab),
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

