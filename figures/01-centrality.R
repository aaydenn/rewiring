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

blood <- tmp2 |> filter(theta == c("Blood.AL","Blood.CCR","Blood.ICR"))

mfp <- tmp2 |> filter(theta == c("MFP.AL","MFP.CCR","MFP.ICR"))

p <- ggpubr::ggboxplot(
  blood,
  "theta",
  "value",
  fill = "centrality",
  bxp.errorbar = T,
  xlab = "coexp. networks",
  ylab = "measurement (scaled)",
  palette = c("gray2", "gray42", "gray82"),
  ggtheme = theme_bw()
) |> ggsave(
  filename = "figures/theta.centrality.blood.png",
  units = "in",
  width = 8,
  height = 8,
  dpi = 300
)

p2 <- ggpubr::ggboxplot(
  blood,
  "theta",
  "value",
  fill = "centrality",
  bxp.errorbar = T,
  xlab = "coexp. networks",
  ylab = "measurement (scaled)",
  palette = c("salmon2", "royalblue2", "springgreen2"),
  ggtheme = theme_bw()
) |> ggsave(
  filename = "figures/theta.centrality.blood_colored.png",
  units = "in",
  width = 8,
  height = 8,
  dpi = 300
)

p3 <- ggpubr::ggboxplot(
  mfp,
  "theta",
  "value",
  fill = "centrality",
  bxp.errorbar = T,
  xlab = "coexp. networks",
  ylab = "measurement (scaled)",
  palette = c("salmon2", "royalblue2", "springgreen2"),
  ggtheme = theme_bw()
) |> ggsave(
  filename = "figures/theta.centrality.mfp_colored.png",
  units = "in",
  width = 8,
  height = 8,
  dpi = 300
)

p4 <- ggpubr::ggboxplot(
  mfp,
  "theta",
  "value",
  fill = "centrality",
  bxp.errorbar = T,
  xlab = "coexp. networks",
  ylab = "measurement (scaled)",
  palette = c("gray2", "gray42", "gray82"),
  ggtheme = theme_bw()
) |> ggsave(
  filename = "figures/theta.centrality.mfp.png",
  units = "in",
  width = 8,
  height = 8,
  dpi = 300
)
