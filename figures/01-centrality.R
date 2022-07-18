library(ggplot2)
library(ggpubr)

load("result/theta_delta.RData")

#' plotting centrality for theta

tmp <- lapply(1:length(theta.centrality), function(x)
  theta.centrality[[x]] |> mutate(
    theta = names(theta.centrality[x]),
    degree = scale(degree),
    knn = scale(knn)
  ) |> select(!hub))

tmp2 <- do.call(rbind, tmp) |> 
  reshape2::melt(variable = "centrality")

blood <- tmp2 |> filter(theta == c("Blood.AL","Blood.CCR","Blood.ICR")) |> mutate(theta=str_remove(theta,"Blood."))

mfp <- tmp2 |> filter(theta == c("MFP.AL","MFP.CCR","MFP.ICR")) |> mutate(theta=str_remove(theta,"MFP."))

p <- ggpubr::ggboxplot(
  blood,
  "theta",
  "value",
  fill = "centrality",
  bxp.errorbar = T,
  xlab = "",
  ylab = "measurement (blood)",
  palette = c("gray2", "gray42", "gray82"),
  ggtheme = theme_bw()
) + theme(text=element_text(size = 18))

ggsave(p,
  filename = "figures/theta.centrality.blood.png",
  units = "in",
  width = 8,
  height = 5,
  dpi = 300
)

p2 <- ggpubr::ggboxplot(
  blood,
  "theta",
  "value",
  fill = "centrality",
  bxp.errorbar = T,
  xlab = "",
  ylab = "measurement (blood)",
  palette = c("salmon2", "royalblue2", "springgreen2"),
  ggtheme = theme_bw()
) + 
  theme(text=element_text(size = 18))

ggsave(p2,
  filename = "figures/theta.centrality.blood_colored.png",
  units = "in",
  width = 8,
  height = 5,
  dpi = 300
)

p3 <- ggpubr::ggboxplot(
  mfp,
  "theta",
  "value",
  fill = "centrality",
  bxp.errorbar = T,
  xlab = "",
  ylab = "measurement (MFP)",
  palette = c("salmon2", "royalblue2", "springgreen2"),
  ggtheme = theme_bw()
) + 
  theme(text=element_text(size = 18))

ggsave(p3,
  filename = "figures/theta.centrality.mfp_colored.png",
  units = "in",
  width = 8,
  height = 5,
  dpi = 300
)

p4 <- ggpubr::ggboxplot(
  mfp,
  "theta",
  "value",
  fill = "centrality",
  bxp.errorbar = T,
  xlab = "",
  ylab = "measurement (MFP)",
  palette = c("gray2", "gray42", "gray82"),
  ggtheme = theme_bw()
) + 
  theme(text=element_text(size = 18))

ggsave(p4,
  filename = "figures/theta.centrality.mfp.png",
  units = "in",
  width = 8,
  height = 5,
  dpi = 300
)



#' plotting centrality for delta

tmp <- lapply(1:length(delta.centrality), function(x)
  delta.centrality[[x]] |> mutate(
    delta = names(delta.centrality[x]),
    degree = scale(degree),
    knn = scale(knn)
  ) |> select(!hub))

tmp2 <- do.call(rbind, tmp) |> 
  reshape2::melt(variable = "centrality")

blood <- tmp2 |> filter(delta == c("D.Blood.CCR.AL","D.Blood.ICR.AL"))

mfp <- tmp2 |> filter(delta == c("D.MFP.CCR.AL","D.MFP.ICR.AL"))


p5 <- ggpubr::ggboxplot(
  blood,
  "delta",
  "value",
  fill = "centrality",
  bxp.errorbar = T,
  xlab = "",
  ylab = "measurement (blood)",
  palette = c("salmon2", "royalblue2", "springgreen2"),
  ggtheme = theme_bw()
) + 
  theme(text=element_text(size = 18))

ggsave(p5,
  filename = "figures/delta.centrality.blood_colored.png",
  units = "in",
  width = 8,
  height = 5,
  dpi = 300
)


p6 <- ggpubr::ggboxplot(
  mfp,
  "delta",
  "value",
  fill = "centrality",
  bxp.errorbar = T,
  xlab = "",
  ylab = "measurement (MFP)",
  palette = c("salmon2", "royalblue2", "springgreen2"),
  ggtheme = theme_bw()
) + 
  theme(text=element_text(size = 18))

ggsave(p6,
  filename = "figures/delta.centrality.mfp_colored.png",
  units = "in",
  width = 8,
  height = 5,
  dpi = 300
)
