#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library("optparse")))
suppressWarnings(suppressMessages(library("limma")))
suppressWarnings(suppressMessages(library("oligo")))
suppressWarnings(suppressMessages(library("dplyr")))
suppressWarnings(suppressMessages(library("readr")))

option_list <- list(
    make_option(
        c("-s", "--sdrf"), type = "character", 
        default = NULL, metavar = "/path/to/sdrf_file/",
        help = "Tab seperated file contains sample and data relations"),
    make_option(
        c("-a", "--annotation"), type="character", 
        default = NULL, metavar = "/path/to/annotation_file/",
        help="Annotation file for microarray")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# read the experiment description 
sdrf <- read_delim(
    opt$sdrf, 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE)

sdrf$`array path` <- paste("Data", sdrf$`organism part`, 
    sdrf$`array data file`, sep = "/")

sdrf |>
  select(`array path`) |>
  read.celfiles(pkgname = "pd.mirna.4.1") -> raw_data

sdrf |> 
    select(diet, `age (week)`, `organism part`) |> 
    mutate(diet = case_when(diet == "ICRR" ~ "ICR", 
        diet == "ICRRF" ~ "ICR", 
        TRUE ~ diet)) |> 
    mutate(diet = as.factor(diet)) |> 
    mutate(diet = relevel(diet, 
        ref="BASELINE")) -> features

colnames(features) <- c("diet", "age", "tissue")

pData(raw_data) <- cbind(pData(raw_data), features)

cat("reading data finished!\n")

# normalize data
eset <- rma(raw_data)
norm_data <- exprs(eset)

cat("normalization finished!\n")

# filter only mouse
read_csv(opt$annotation, comment = "#", show_col_types = FALSE) |> 
    filter(
    `Species Scientific Name` == "Mus musculus", 
    # only mature miRNAs
    str_starts(Accession, "MIMAT")
    ) -> mouse_annot

mouse_data <- norm_data[mouse_annot$`Probe Set Name`, ]
rownames(mouse.Data) <- mouse_annot$`Transcript ID(Array Design)`

save(mouse_data, features, file = "mouse_data_features.rda")

if(nrow(mouse_data)>0) {
    cat("filtering finished!\n")
} else {
    cat("filtering gives zero value!\n")
    break
}

# PCA
raw_data |> exprs() |> 
    t() |> log1p() |> prcomp() -> raw_pca

norm_data |> 
    t() |> log1p() |> prcomp() -> norm_pca

raw_var <- round(100 * raw_pca$sdev ^ 2 / sum(raw_pca$sdev ^ 2), 1)
norm_var <- round(100 * norm_pca$sdev ^ 2 / sum(norm_pca$sdev ^ 2), 1)

# data frame for ggplot aesthetics
raw_data_frame <- data.frame(
    PC1 = raw_pca$x[,1], 
    PC2 = raw_pca$x[,2], 
    Tissue = pData(raw_data)$tissue,
    Diet = pData(raw_data)$diet)

norm_data_frame <- data.frame(
    PC1 = norm_pca$x[,1], 
    PC2 = norm_pca$x[,2], 
    Tissue = features$tissue,
    Diet = features$diet)

p_raw <- 
    ggplot(raw_data_frame, aes(PC1, PC2)) + 
    geom_point(aes(colour = Tissue, shape = Diet)) +
    xlab(paste0("PC1, VarExp: ", raw_var[1], "%")) +
    ylab(paste0("PC2, VarExp: ", raw_var[2], "%")) +
    scale_color_manual(
        values = c("#141542", '#a7a523', '#809e52')
    ) + 
    ggtitle("PCA plot of the log-transformed raw expression data") + 
    theme_bw()

ggsave("PCAplotRaw.png", p_raw, units = "px", 
    width = 1725, height = 1112, dpi = 300)

p_norm <- 
    ggplot(norm_data_frame, aes(PC1, PC2)) + 
    geom_point(aes(colour = Tissue, shape = Diet)) +
    xlab(paste0("PC1, VarExp: ", norm_var[1], "%")) +
    ylab(paste0("PC2, VarExp: ", norm_var[2], "%")) +
    scale_color_manual(
        values = c("#141542", '#a7a523', '#809e52')
    ) + 
    ggtitle("PCA plot of the log-transformed raw expression data") + 
    theme_bw()

ggsave("PCAplotNorm.png", p_norm, units = "px", 
    width = 1725, height = 1112, dpi = 300)

cat("PCA plots created!\n")