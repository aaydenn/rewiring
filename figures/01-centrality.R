library(ggplot2)
library(ggpubr)

load("result/theta_delta.RData")



#' plotting centrality for theta

tmp <- lapply(1:length(theta.centrality), function(x)
  theta.centrality[[x]] |> mutate(
    theta = names(theta.centrality[x]),
    degree = degree,
    knn = knn
  ) |> select(!hub))


tmp2 <- do.call(rbind, tmp) |> 
  reshape2::melt(variable = "centrality")




blood <- tmp2 |> 
  filter(theta == c("Blood.AL","Blood.CCR","Blood.ICR")) |> 
  mutate(theta=str_remove(theta,"Blood."))



p <- ggplot(data = blood, aes(x = theta, y = value, fill = centrality)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("salmon2", "royalblue2", "springgreen2")) + 
  facet_wrap(~centrality,scales = "free_y") + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.title.x = element_blank(), 
        legend.position = "none") + 
  ylab("centrality indices")


ggsave(p, filename = "figures/theta.centrality.blood.png",
  units = "in", width = 8, height = 5, dpi = 300)





mfp <- tmp2 |> 
  filter(theta == c("MFP.AL","MFP.CCR","MFP.ICR")) |> 
  mutate(theta=str_remove(theta,"MFP."))



p1 <- ggplot(data = mfp, aes(x = theta, y = value, fill = centrality)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("salmon2", "royalblue2", "springgreen2")) + 
  facet_wrap(~centrality, scales = "free_y") + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.title.x = element_blank(), 
        legend.position = "none") + 
  ylab("centrality indices")

ggsave(p1, filename = "figures/theta.centrality.mfp_colored.png",
  units = "in", width = 8, height = 5, dpi = 300)





#' plotting centrality for delta

tmp0 <- lapply(1:length(delta.centrality), function(x)
  delta.centrality[[x]] |> mutate(
    delta = names(delta.centrality[x]),
    degree = degree,
    knn = knn) |> 
    select(!hub))

tmp02 <- do.call(rbind, tmp0) |> 
  reshape2::melt(variable = "centrality")



blood0 <- tmp02 |> filter(delta == c("D.Blood.CCR.AL","D.Blood.ICR.AL"))  |> 
  mutate(delta = str_remove(delta,"D.Blood."))



p2 <- ggplot(data = blood0, aes(x = delta, y = value, fill = centrality)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("salmon2", "royalblue2", "springgreen2")) + 
  facet_grid(~centrality) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.title.x = element_blank(), 
        legend.position = "none") + 
  ylab("centrality indices")


ggsave(p2, filename = "figures/delta.centrality.blood_colored.png",
  units = "in", width = 8, height = 5, dpi = 300)





mfp0 <- tmp02 |> filter(delta == c("D.MFP.CCR.AL","D.MFP.ICR.AL")) |> 
  mutate(delta = str_remove(delta,"D.MFP."))



p3 <- ggplot(data = mfp0, aes(x = delta, y = value, fill = centrality)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("salmon2", "royalblue2", "springgreen2")) + 
  facet_grid(~centrality) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.title.x = element_blank(), 
        legend.position = "none") + 
  ylab("centrality indices")

ggsave(p3, filename = "figures/delta.centrality.mfp_colored.png",
  units = "in", width = 8, height = 5, dpi = 300)
