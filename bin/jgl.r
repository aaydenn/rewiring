#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library("DiffGraph")))
suppressWarnings(suppressMessages(library("optparse")))
suppressWarnings(suppressMessages(library("dplyr")))
suppressWarnings(suppressMessages(library("igraph")))
suppressWarnings(suppressMessages(library("randomcoloR")))


option_list <- list(
    make_option(
        c("-n", "--normalized-data"), type = "character", 
        default = NULL, metavar = "/path/to/rdata/",
        help = "R data matrix object that contains normalized microarray data"),
    make_option(
        c("-f", "--features"), type = "character", 
        default = NULL, metavar = "/path/to/rdata/",
        help = "R data matrix object that contains features of the samples"),
    make_option(
        c("-m", "--mirna"), type = "character", 
        default = NULL, metavar = "/path/to/rdata/",
        help = "R data vector object(s) that contains significant miRNAs")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

load(opt$n)
load(opt$f)
load(opt$m)

mouse_data |>
  as.data.frame() |> 
  filter(rownames(mouse_data) %in% c(de_blood, de_mfp)) |> 
  t() -> arrays

colnames(arrays) <- str_remove(colnames(arrays), "mmu-")

if(nrow(arrays)>0) {
    cat("arrays are selected\n")
} else {
    cat("empty arrays are selected!\n")
}

# function for centrality measures for given network
centrality_table <- function(net) { 
    data.frame(
        degree = strength(net), 
        closeness = closeness(net) |> 
            signif(digits=3),
        betweenness = betweenness(net) |> 
            signif(digits=3),
        eigen = eigen_centrality(net)$vector |> 
            signif(digits=3),
        knn = knn(net, V(net))$knn |> 
            signif(digits=3)) |> 
        rownames_to_column(var = "miRNA")
}

# function for network comparison of two diets for given tissue
fgl_for_diet_tissue <- function(diet1, diet2, tissue) {
    arrays_diet1 <- 
        arrays[features$tissue == tissue & features$diet %in% c(diet1, "BASELINE"),]
    arrays_diet2 <- 
        arrays[features$tissue == tissue & features$diet %in% c(diet2, "BASELINE"),]
    # combinations for Differential Co-expression (DCE) graphs
    d_diets <- 
        setNames(list(arrays_diet1, arrays_diet2) , 
            c(diet1, diet2)) 
    d_fgl <- 
        FGL(d_diets, 
            lambda1 = 0.05, lambda2 = 0.25, 
            covType = "spearman")
    d_centrality <- 
        centrality_table(d_fgl$Delta.graph.connected)
    return(
        list(fgl = d_fgl, centrality = d_centrality)
    )
}

diets <- c(CCR = "CCR", ICR = "ICR")
tissues <- c(Blood = "Blood", MFP = "MFP")

d_fgl <- lapply(
    diets, 
        function(y) lapply(
            tissues, 
                function(x) fgl_for_diet_tissue(y, "AL", x))
)

save(d_fgl, file = "d_fgl.rda")

cat("networks constructed!\n")

for(diet in diets) {
    for(tissue in tissues) {
        write.table(x = d_fgl[[diet]][[tissue]]$centrality, 
            file = paste(diet, "AL" , tissue, "centrality.tsv", sep = "_"))
  }
}

cat("centralities calculated!\n")

find_communities <- function(diet, tissue) {
    d_net <- d_fgl[[diet]][[tissue]]$fgl$Delta.graph.connected
    cluster_leiden(d_net, resolution_parameter = 0.1)
}

communities <- lapply(
    diets, 
        function(diet) lapply(
            tissues, 
                function(tissue) find_communities(diet, tissue))
)

cat("communities found!\n")

draw_net_diet_tissue <- function(diet, tissue) {
    d_net <- 
        d_fgl[[diet]][[tissue]]$fgl$Delta.graph.connected
    d_com <- 
        communities[[diet]][[tissue]]
    communityColors <- 
        distinctColorPalette(max(d_com$membership))[membership(d_com)]

    png(paste(diet, "AL" , tissue, "net.png", sep = "_"), 10, 10, units = "in", res = 320)
    plot(
        d_com, d_net, mark.groups = NULL, vertex.label.cex = 0.5,
        vertex.size = degree(d_net) * 1.5,
        col = communityColors, layout = layout_with_kk,
        main = paste("DCE network of", diet, "in", tissue)
    )
    dev.off()
}

for(diet in diets) {
  for(tissue in tissues) {
    draw_net_diet_tissue(diet, tissue)
  }
}

cat("networks printed!\n")