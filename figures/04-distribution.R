library(igraph)

load("Data/fgl.Rdata") # first level mirna-mirna (GLASSO) ####
load("Data/targets.Rdata") # second level mirna-gene (mirnatap) ####

#' ---------------------------------
# degree distribution for each delta

delta.degree.dist <- lapply(delta.graphs, function(x) {
  x |> 
    degree() |> 
    table() |> 
    as.data.frame() |> 
    mutate(Var1 = as.numeric(Var1))
})

blood <- rbind(
  data.frame(delta.degree.dist$D.Blood.CCR.AL, `Network ID` = "D.Blood.CCR.AL"),
  data.frame(delta.degree.dist$D.Blood.ICR.AL, `Network ID` = "D.Blood.ICR.AL")
)

p1 <- ggplot(blood, aes(Var1, log10(Freq), colour=Network.ID)) + 
  geom_point()+
  geom_smooth(se=FALSE)+
  scale_color_manual(values=c("salmon2", "royalblue2"))+
  geom_ribbon(stat = "smooth",
              method = "loess",
              se = TRUE,
              alpha = 0.1, # or, use fill = NA
              colour = "black",
              linetype = "dotted")+
  xlab("Degree(k)") + 
  ylab("log(Frequency)") +
  theme_bw()+
  theme(legend.position = c(0.9, 0.9), legend.justification = c(1,1), legend.background = element_blank())+
  theme(text=element_text(family="Times", size = 18));p1

mfp <- rbind(
  data.frame(delta.degree.dist$D.MFP.CCR.AL, `Network ID` = "D.MFP.CCR.AL"),
  data.frame(delta.degree.dist$D.MFP.ICR.AL, `Network ID` = "D.MFP.ICR.AL")
)

p2 <- ggplot(mfp, aes(Var1, log10(Freq), colour=Network.ID)) + 
  geom_point()+
  geom_smooth(se=FALSE)+
  scale_color_manual(values=c("salmon2", "royalblue2"))+
  geom_ribbon(stat = "smooth",
              method = "loess",
              se = TRUE,
              alpha = 0.1, # or, use fill = NA
              colour = "black",
              linetype = "dotted")+
  xlab("Degree(k)") + 
  ylab("log(Frequency)") +
  theme_bw()+
  theme(legend.position = c(0.9, 0.9), legend.justification = c(1,1), legend.background = element_blank())+
  theme(text=element_text(family="Times", size = 18));p2

ggsave("figures/delta.degree.dist.blood.png", p1, width = 8, height = 8, units = "in", dpi = 300)
ggsave("figures/delta.degree.dist.mfp.png",   p2, width = 8, height = 8, units = "in", dpi = 300)
