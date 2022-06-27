# Principle Component Analysis

library(Biobase)

FULL_DATA <- "Data/full.Data.RData"

#####
# use full data (unfiltered) for pca
#####
load(FULL_DATA)

windowsFonts(Times = windowsFont("Times New Roman"))

# Before normalizing ####
## expression data is commonly analyzed
## on a logarithmic scale (log2 of full.Data)
exp_raw0 <- log2(Biobase::exprs(full.Data))

pData(full.Data) <- cbind(pData(full.Data), features)
rownames(features) <- colnames(exp_raw0)

## ignore brain data
keep <- features[features$tissue != "Brain", ]
exp_raw <- exp_raw0[, rownames(keep)]

p_data0 <- pData(full.Data)[pData(full.Data)$tissue != "Brain", ]

# TODO: fix diet names earlier
p_data <- mutate(p_data0,
    diet = ifelse(diet == "ICRR" | diet == "ICRRF", "ICR", diet))

## calculate principal components
PCAraw <- prcomp(t(exp_raw), scale. = FALSE)
percentVarRaw <- round(100 * PCAraw$sdev^2 / sum(PCAraw$sdev^2), 1)
sd_ratio_raw <- sqrt(percentVarRaw[2] / percentVarRaw[1])

# data frame for ggplot aesthetics
dataGGraw <- data.frame(PC1 = PCAraw$x[,1],
                        PC2 = PCAraw$x[,2],
                        Tissue = p_data$tissue,
                        Diet = p_data$diet)
# make the 2d plot of components
p <- ggplot(dataGGraw, aes(PC1, PC2, size = 10)) +
  geom_point(aes(colour = Tissue,
                 shape  = Diet)) +
  xlab(paste0("PC1, VarExp: ", percentVarRaw[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVarRaw[2], "%")) +
  scale_color_manual(values = c("salmon2", 'royalblue2')) +
  # coord_fixed(ratio = sd_ratio_raw*2) + 
  ggtitle("PCA plot of the log-transformed raw expression data") + 
  theme_bw() +
  theme(text=element_text(family="Times",size = 18))
p# save the plot
ggsave("PCAplotRaw.png", p, units = "in", 
       width = 10, height = 10, dpi = 300)
# After normalizing ####
exp_norm0 <- Biobase::exprs(eset)
exp_norm <- exp_norm0[,rownames(keep)]
PCA_norm <- prcomp(t(exp_norm), scale. = FALSE)
percentVarNorm <- round(100*PCA_norm$sdev^2/sum(PCA_norm$sdev^2),1)
sd_ratio_norm <- sqrt(percentVarNorm[2] / percentVarNorm[1])
# data frame for ggplot aesthetics
dataGGnorm <- data.frame(PC1 = PCA_norm$x[,1], 
                         PC2 = PCA_norm$x[,2], 
                         Tissue = p_data$tissue,
                         Diet   = p_data$diet)
# make the 2d plot of components
p0 <- ggplot(dataGGnorm, aes(PC1, PC2, colour = Tissue,
                            shape  = Diet)) + 
  geom_point(size=2) + 
  xlab(paste0("PC1, VarExp: ", percentVarNorm[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVarNorm[2], "%")) +
  scale_color_manual(values = c("salmon2", "royalblue2")) + theme_bw()
  #   +
  # theme(legend.position = c(0.99, 0.01), legend.justification = c(1,0),
  #       legend.background = element_blank())
  # coord_fixed(ratio = sd_ratio_norm) + 
  # theme(legend.position = "bottom") + 

# +
# theme(text=element_text(family="Times", size = 18))
p0
ggsave("PCAplotNorm.png", p0, units = "in", 
       width = 5, height = 4, dpi = 300)
#####
# Use mouse.data for pca
#####
load("Data/mouse-data-features.RData")
## make phenotype data to assign conditions to arrays
fdata <- features
rownames(fdata) <- colnames(mouse.Data)
keep <- fdata[fdata$tissue!="Brain",]
rownames(keep)
mouse.Data <- mouse.Data[,rownames(keep)]
fdata <- fdata[fdata$tissue != "Brain",]
fdata <- mutate(fdata, diet = ifelse(diet == "ICRR" | diet == "ICRRF", "ICR", diet))
# pData(full.Data) <- cbind(pData(full.Data),features)
## calculate principal components
PCA.filtered <- prcomp(t(mouse.Data), scale. = FALSE)
percentVarFiltered <- round(100*PCA.filtered$sdev^2/sum(PCA.filtered$sdev^2),1)
# sd_ratio_raw <- sqrt(percentVarRaw[2] / percentVarRaw[1])
# data frame for ggplot aesthetics
dataGGfiltered <- data.frame(PC1 = PCA.filtered$x[,1], 
                        PC2 = PCA.filtered$x[,2], 
                        Tissue = fdata$tissue,
                        Diet = fdata$diet)
# make the 2d plot of components
p <- ggplot(dataGGfiltered, aes(PC1, PC2)) + 
  geom_point(aes(colour = Tissue,
                 shape  = Diet)) + 
  xlab(paste0("PC1, VarExp: ", percentVarFiltered[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVarFiltered[2], "%")) +
  scale_color_manual(values = c("#141542", '#a7a523', '#809e52')) +
  # coord_fixed(ratio = sd_ratio_norm) + 
  # theme(legend.position = "bottom") + 
  ggtitle("PCA plot of the log-transformed normalized expression data");p
ggsave("PCAplotfiltered.png", p, units = "px", 
       width = 1725, height = 1112, dpi = 300)
