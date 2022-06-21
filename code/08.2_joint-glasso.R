library(DiffGraph)

load("Data/mouse-data-features.RData")
load("Data/DE_miRNAs.RData")

# modify diet column of features
features$diet <- ifelse(features$diet == "ICRR",  "ICR", features$diet)
features$diet <- ifelse(features$diet == "ICRRF", "ICR", features$diet)

mouse.Data |> # filter miRNAs according to DE miRNAs not in Brain
  as.data.frame() |> 
  filter(rownames(mouse.Data) %in% c(DE.blood, DE.mfp)) |> t() -> arrays # filter(rownames(mouse.Data) %in% gnames.unique) |> t() -> arrays

colnames(arrays) <- str_remove(colnames(arrays), "mmu-")
colnames(features) <- c("diet", "age", "tissue", "array.design")

# DIET COMPARISONS ####
centrality.table <- function(net) { # function for centrality measures for given network
  data.frame(degree = strength(net), closeness = closeness(net) |> signif(digits=3),
             betweenness = betweenness(net) |> signif(digits=3),
             eigen = eigen_centrality(net)$vector |> signif(digits=3),
             knn = knn(net, V(net))$knn |> signif(digits=3)) |> 
    rownames_to_column(var = "miRNA")
}


fgl_for_diet_tissue <- function(diet1, diet2, tissue) { # function for network comparison of two diets for given tissue
  arrays.diet1 <- arrays[features$tissue == tissue & features$diet %in% c(diet1, "BASELINE"),]
  arrays.diet2 <- arrays[features$tissue == tissue & features$diet %in% c(diet2, "BASELINE"),]
  d.diets <- setNames(list(arrays.diet1, arrays.diet2) , c(diet1, diet2)) # combinations for Differential Co-expression (DCE) graphs
  d.fgl <- FGL(d.diets, lambda1 = 0.05, lambda2 = 0.25, covType = "spearman")
  d.centrality <- centrality.table(d.fgl$Delta.graph.connected)
  return(list(arrays.diet1 = arrays.diet1, arrays.diet2 = arrays.diet2, 
              d.fgl = d.fgl, d.centrality = d.centrality))
}
# fgl_result <- list()
# k <- 1
# for (i in c("Blood", "Brain", "MFP")) {
#   for (j in c("CCR","ICR")) {
#     fgl_result[[k]] <- fgl_for_diet_tissue(j, "AL", i)
#     k <- k + 1
#   }
# }

# iteratively make fgl to each tissue and diet
d.ccr.al.fgl <- setNames(lapply(c("Blood", "Brain", "MFP"), function(x) fgl_for_diet_tissue("CCR", "AL", x)), c("Blood", "Brain", "MFP"))
d.icr.al.fgl <- setNames(lapply(c("Blood", "Brain", "MFP"), function(x) fgl_for_diet_tissue("ICR", "AL", x)), c("Blood", "Brain", "MFP"))

d.fgl <- setNames(lapply(c("CCR", "ICR"), function(y)
  setNames(lapply(c("Blood", "Brain", "MFP"), function(x)
      fgl_for_diet_tissue(y, "AL", x)),
    c("Blood", "Brain", "MFP"))), c("CCR", "ICR"))

for(diet1 in c("CCR","ICR")) {
  for(tissue in c("Blood", "Brain", "MFP")) {
    write.table(x = d.fgl[[diet1]][[tissue]]$d.centrality,
                file = paste("Results/d", diet1, "AL" , tissue, "centrality.tsv", sep = "."))
  }
}

save(d.ccr.al.fgl, d.icr.al.fgl, d.fgl, file = "Data/d.fgl.RData")

thetagraphs0 <- list(
  al.blood = d.fgl$CCR$Blood$d.fgl$Theta.graph.connected[[2]],
  ccr.blood = d.fgl$CCR$Blood$d.fgl$Theta.graph.connected[[1]],
  icr.blood = d.fgl$ICR$Blood$d.fgl$Theta.graph.connected[[1]],
  al.mfp = d.fgl$CCR$MFP$d.fgl$Theta.graph.connected[[2]],
  ccr.mfp = d.fgl$CCR$MFP$d.fgl$Theta.graph.connected[[1]],
  icr.mfp = d.fgl$ICR$MFP$d.fgl$Theta.graph.connected[[1]]
)

theta.centrality1 <- lapply(thetagraphs0, centrality.table)

theta.centrality <- lapply(theta.centrality1, function(x) x |> mutate(hub = ifelse(knn<degree, 1,0)))
# openxlsx::write.xlsx(theta.centrality,"Results/theta.centrality.xlsx")


