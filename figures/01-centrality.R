library(ggplot2)
library(ggpubr)

load("result/fgl.RData")

#' plotting centrality

tmp <- lapply(1:length(theta.centrality), function(x)
  theta.centrality[[x]] |> mutate(
    theta = names(theta.centrality[x]),
    degree = scale(degree),
    knn = scale(knn)
  ) |> select(!hub))

tmp2 <- do.call(rbind, tmp) |> 
  reshape2::melt(variable = "centrality")

p <- ggpubr::ggboxplot(
  tmp2,
  "theta",
  "value",
  fill = "centrality",
  bxp.errorbar = T,
  xlab = "coexp. networks",
  ylab = "measurement (scaled)",
  palette = c("gray2", "gray42", "gray82"),
  ggtheme = theme_bw()
)

ggsave("figures/theta.centrality.png", p, units = "in", 
       width = 5, height = 4, dpi = 300)