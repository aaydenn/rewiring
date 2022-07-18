library(stringr)
library(ggpubr)

load("result/mouse-data-features.RData")

colnames(mouse.Data) <- toupper(colnames(mouse.Data))
colnames(mouse.Data) <- gsub("ICRR|ICRRF", "ICR", colnames(mouse.Data))

as.data.frame(mouse.Data) |> 
  # T95: blood and mfp same diet, age
  ggscatter(x = "AL_T95_BLOOD_W50", 
            y = "AL_T95_MFP_W50",
            add = "reg.line", conf.int = TRUE, 
            cor.coef = F, cor.method = "pearson",
            xlab = "Blood of AL fed mice", 
            ylab = "MFP of AL fed mice",
            ggtheme = theme_bw(),
            color	= alpha("black", 0.5)) + 
  theme(text=element_text(size = 18)) -> p1

ggsave(p1, file = "figures/same_animal_T95_AL.png", width = 8, height = 8, units = "in", dpi = 300)

as.data.frame(mouse.Data) |> 
  # T95 and T106: different diet, same tissue, age
  ggscatter(x = "AL_T95_MFP_W50",
            y = "CCR_T106_MFP_W50",
            add = "reg.line", conf.int = TRUE, 
            cor.coef = F, cor.method = "pearson",
            xlab = "MFP of AL fed mice", 
            ylab = "MFP of CCR fed mice",
            ggtheme = theme_bw(),
            color	= alpha("black", 0.5)
  ) + theme(text=element_text(size = 18)) -> p2

ggsave(p2, file = "figures/same_diet_T95_T106_mfp.png", width = 8, height = 8, units = "in", dpi = 300)

as.data.frame(mouse.Data) |>
  # T95 and T106: different diet, same tissue, age
  ggscatter(x = "AL_T95_BLOOD_W50",
            y = "CCR_T106_BLOOD_W50",
            add = "reg.line", conf.int = TRUE, 
            cor.coef = F, cor.method = "pearson",
            xlab = "Blood of AL fed mice", 
            ylab = "Blood of CCR fed mice",
            ggtheme = theme_bw(),
            color	= alpha("black", 0.5)) + theme(text=element_text(size = 18)) -> p3

ggsave(p3,file = "figures/same_diet_T95_T106_blood.png", width = 8, height = 8, units = "in", dpi = 300)

