library(stringr)
library(ggpubr)
library(Biobase)
library(patchwork)

load("result/mouse-data-features.RData")
load("G:/Drive'ım/CR-miRNA/Data/full.Data.RData")

colnames(mouse.Data) <- toupper(colnames(mouse.Data))
colnames(mouse.Data) <- gsub("ICRR|ICRRF", "ICR", colnames(mouse.Data))

size <- 12

# Correlations

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
  theme(text=element_text(size = size)) -> p1

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
            color	= alpha("black", 0.5)) + 
  theme(text=element_text(size = size)) -> p2

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
            color	= alpha("black", 0.5)) + 
  theme(text=element_text(size = size)) -> p3

ggsave(p3,file = "figures/same_diet_T95_T106_blood.png", width = 8, height = 8, units = "in", dpi = 300)


p4 <- p2 + p3

ggsave("figures/same_diet.png",p4,width = 8, height = 4, units = "in", dpi = 300)


# PCA
# Raw expression


pData(full.Data) <- cbind(pData(full.Data), features)



## ignore brain data



p_data0 <- pData(full.Data)[pData(full.Data)$tissue != "Brain", ]

p_data <- mutate(p_data0, 
  diet = ifelse(diet == "ICRR" | diet == "ICRRF", "ICR", diet))


# Before normalizing

exp_raw0 <- log2(Biobase::exprs(full.Data))

rownames(features) <- colnames(exp_raw0)

keep <- features[features$tissue != "Brain", ]

exp_raw <- exp_raw0[, rownames(keep)]



## calculate principal components

PCAraw <- prcomp(t(exp_raw), scale. = FALSE)

percentVarRaw <- round(100 * PCAraw$sdev^2 / sum(PCAraw$sdev^2), 1)

# sd_ratio_raw <- sqrt(percentVarRaw[2] / percentVarRaw[1])



# data frame for ggplot aesthetics

dataGGraw <- data.frame(PC1 = PCAraw$x[,1], PC2 = PCAraw$x[,2],
                        Tissue = p_data$tissue, Diet = p_data$diet)


p01 <- ggplot(dataGGraw, aes(PC1, PC2)) +
  geom_point(aes(colour = Tissue,
                 shape  = Diet),size = 5) +
  xlab(paste0("PC1, VarExp: ", percentVarRaw[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVarRaw[2], "%")) +
  scale_color_manual(values = c("salmon2", 'royalblue2')) +
  # coord_fixed(ratio = sd_ratio_raw*2) + 
  # ggtitle("PCA plot of the log-transformed raw expression data") + 
  theme_bw() +
  theme(# legend properties
    legend.position = c(0.99, 0.01), 
    legend.justification = c(1,0), 
    legend.background = element_rect(fill = alpha("grey", 0.2)),
    legend.key.height = unit(0.1,"in"),
    legend.key.size = unit(0.1,"in"),
    # every text element
    text = element_text(size = size))


ggsave("figures/pca_raw.png", p01, units = "in", width = 8, height = 8, dpi = 300)


# PCA
# Norm. expression


exp_norm0 <- Biobase::exprs(eset)

rownames(features) <- colnames(exp_norm0)

keep <- features[features$tissue != "Brain", ]

exp_norm <- exp_norm0[,rownames(keep)]



PCA_norm <- prcomp(t(exp_norm), scale. = FALSE)

percentVarNorm <- round(100*PCA_norm$sdev^2/sum(PCA_norm$sdev^2),1)
# sd_ratio_norm <- sqrt(percentVarNorm[2] / percentVarNorm[1])


# data frame for ggplot aesthetics

dataGGnorm <- data.frame(PC1 = PCA_norm$x[,1], PC2 = PCA_norm$x[,2], 
                         Tissue = p_data$tissue, Diet   = p_data$diet)


p02 <- ggplot(dataGGnorm, aes(PC1, PC2)) + 
  geom_point(aes(colour = Tissue,
             shape  = Diet), size = 3) + 
  xlab(paste0("PC1, VarExp: ", percentVarNorm[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVarNorm[2], "%")) +
  scale_color_manual(values = c("salmon2", "royalblue2")) + 
  # coord_fixed(ratio = sd_ratio_raw*2) + 
  # ggtitle("PCA plot of the log-transformed norm. expression data") + 
  theme_bw() +
  theme(# legend properties
    legend.position = c(0.99, 0.01), 
    legend.justification = c(1,0), 
    legend.background = element_rect(fill = alpha("grey", 0.2)),
    legend.key.height = unit(0.1,"in"),
    legend.key.size = unit(0.1,"in"),
    # every text element
    text = element_text(size = size))

ggsave("figures/pca_norm.png", p02, units = "in", width = 8, height = 8, dpi = 300)


# for manuscript

((p02 + p1 + p2 + p3) + plot_annotation(tag_levels = 'A')) |>
  ggsave(filename = "similarity_manuscript.svg", width = 8, height = 8, dpi = 320)
