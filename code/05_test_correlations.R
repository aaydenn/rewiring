library(stringr)
library(ggpubr)
source("Code/chart.Correlation.R")
#library(limma)

load("Data/mouse-data-features.RData")

colnames(mouse.Data) <- toupper(colnames(mouse.Data))
colnames(mouse.Data) <- gsub("ICRR|ICRRF", "ICR", colnames(mouse.Data))

to_show <- c("AL_T95_BRAIN_W50", "CCR_T106_BRAIN_W50", "ICR_T116_BRAIN_W80", 
             "AL_T95_BLOOD_W50", "CCR_T106_BLOOD_W50", "ICR_T116_BLOOD_W80",
             "AL_T95_MFP_W50", "CCR_T106_MFP_W50", "ICR_T116_MFP_W80")

plot(as.data.frame(mouse.Data[,to_show]), pch=".")
plot(as.data.frame(mouse.Data[,to_show[c(1,4,7)]]), pch=".", main="Same animal")
plot(as.data.frame(mouse.Data[,to_show[1:3]]), pch=".", main="Brain")
plot(as.data.frame(mouse.Data[,to_show[4:6]]), pch=".", main="Blood")
plot(as.data.frame(mouse.Data[,to_show[7:9]]), pch=".", main="MFP")
cor(mouse.Data[,to_show]) |> round(3)
cor(mouse.Data[,to_show]) |>
  heatmap(Rowv = NA, Colv = NA, scale = "none", revC = TRUE)

# visualization ####
png( 
  filename = "Results/same_animal_cor.png",
  res = 300,
  units = "in",
  width = 10,
  height = 10
)
chart.Correlation( # same animal all tissues
  as.data.frame(mouse.Data[, to_show[c(1, 4, 7)]]),
  histogram = T,
  method = "pearson",
  main = "Same Animal"
)
dev.off()

png( 
  filename = "Results/same_tissue_blood_cor.png",
  res = 300,
  units = "in",
  width = 10,
  height = 10
)
chart.Correlation( # different animals same tissue (blood) 
  as.data.frame(mouse.Data[, to_show[4:6]]),
  histogram = T,
  method = "pearson",
  main = "Same Animal"
)
dev.off()

png(
  filename = "Results/same_tissue_brain_cor.png",
  res = 300,
  units = "in",
  width = 10,
  height = 10
)
chart.Correlation(
  as.data.frame(mouse.Data[, to_show[1:3]]),
  histogram = T,
  method = "pearson",
  main = "Same Animal"
)
dev.off()

png(
  filename = "Results/same_tissue_mfp_cor.png",
  res = 300,
  units = "in",
  width = 10,
  height = 10
)
chart.Correlation(
  as.data.frame(mouse.Data[, to_show[7:9]]),
  histogram = T,
  method = "pearson",
  main = "Same Animal"
)
dev.off()

png(
  filename = "Results/cor.png",
  res = 300,
  units = "in",
  width = 10,
  height = 10
)
chart.Correlation(as.data.frame(mouse.Data[, to_show]),
                  histogram = T,
                  method = "pearson")
dev.off()

png(
  filename = "Results/same_animal_mfp_blood_cor.png",
  res = 300,
  units = "in",
  width = 10,
  height = 10
)
chart.Correlation(as.data.frame(mouse.Data[, to_show[c(4,7)]]),
                  histogram = T,
                  method = "pearson")
dev.off()

###############################################################################
# annotation <- features[which(rownames(features) %in% to_show), ]
# pheatmap(
#   cor(mouse.Data[, to_show]),
#   annotation_row = annotation[, c(1, 3)],
#   cluster_rows = F,
#   cluster_cols = F,
#   show_rownames = F,
#   display_numbers = T,
#   filename = "Results/heatmap_cor.png",
#   res = 300,
#   units = "in",
#   width = 10,
#   height = 10
# )
###############################################################################

as.data.frame(mouse.Data) |> 
  # T95: blood and mfp same diet, age
  ggscatter(x = "AL_T95_BLOOD_W50", 
            y = "AL_T95_MFP_W50",
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Blood of AL fed mice", 
            ylab = "MFP of AL fed mice",
            ggtheme = theme_bw(),
            color	= alpha("black", 0.5)) -> p1

ggsave(p1, file = "same_animal_T95_AL.png", width = 6, height = 6, units = "in", dpi = 300)

as.data.frame(mouse.Data) |> 
  # T106: blood and mfp same diet, age
  ggscatter(x = "CCR_T106_BLOOD_W50", 
            y = "CCR_T106_MFP_W50",
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Blood of CCR fed mice (T106)", 
            ylab = "MFP of CCR fed mice (T106)",
            ggtheme = theme_bw(),
            title = "Expression Profiles Between Different Tissues in Same Animal") |> 
  ggpar(font.title = c(18, "black"),
        font.x = c(18,  "black"),
        font.y = c(18,  "black")) |> 
  ggsave(file = "same_animal_T106_AL.png", width = 10, height = 10, units = "in", dpi = 320)

#####

as.data.frame(mouse.Data) |> 
  # T95 and T106: different diet, same tissue, age
  ggscatter(x = "AL_T95_MFP_W50",
            y = "CCR_T106_MFP_W50",
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "MFP of AL fed mice (T95)", 
            ylab = "MFP of CCR fed mice (T106)",
            ggtheme = theme_bw(),
            color	= alpha("black", 0.5)
            # title = "Expression Profiles Between Different Diets in Same Tissue"
            ) -> p2

p2
ggsave(p2, file = "same_diet_T95_T106_mfp.png", width = 10, height = 10, units = "in", dpi = 320)


as.data.frame(mouse.Data) |> 
  # T95 and T116: different diet, same tissue, age
  ggscatter(x = "AL_T95_MFP_W50",
            y = "ICR_T116_MFP_W80",
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "MFP of AL fed mice (T95)", 
            ylab = "MFP of ICR fed mice (T116)",
            ggtheme = theme_bw(),
            title = "Expression Profiles Between Different Diets in Same Tissue") |> 
  ggpar(font.title = c(18, "black"),
        font.x = c(18,  "black"),
        font.y = c(18,  "black")) |> 
  ggsave(file = "same_diet_T95_T116_mfp.png", width = 10, height = 10, units = "in", dpi = 320)


#####

as.data.frame(mouse.Data) |>
  # T95 and T106: different diet, same tissue, age
  ggscatter(x = "AL_T95_BLOOD_W50",
            y = "CCR_T106_BLOOD_W50",
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "BLOOD of AL fed mice (T95)", 
            ylab = "BLOOD of CCR fed mice (T106)",
            ggtheme = theme_bw(),
            color	= alpha("black", 0.5)) -> p3

  ggsave(file = "same_diet_T95_T106_blood.png", width = 10, height = 10, units = "in", dpi = 320)


as.data.frame(mouse.Data) |> 
  # T95 and T116: different diet, same tissue, age
  ggscatter(x = "AL_T95_BLOOD_W50",
            y = "ICR_T116_BLOOD_W80",
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "BLOOD of AL fed mice (T95)", 
            ylab = "BLOOD of ICR fed mice (T116)",
            ggtheme = theme_bw(),
            title = "Expression Profiles Between Different Diets in Same Tissue") |> 
  ggpar(font.title = c(18, "black"),
        font.x = c(18,  "black"),
        font.y = c(18,  "black")) |> 
  ggsave(file = "same_diet_T95_T116_blood.png", width = 10, height = 10, units = "in", dpi = 320)


ggarrange(p0, p1, p2, p3, labels=LETTERS[1:4], ncol = 2, nrow = 2) |>
  ggsave(file="test.png", width = 8, height = 8, units = "in", dpi = 300)


