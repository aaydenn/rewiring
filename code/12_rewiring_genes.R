library(igraph)
library(graphite)
library(pbapply)
library(ggrepel)

# load("Data/d.fgl.RData")
load("Data/targets.Rdata")
# load("Data/miRNA_influence.RData")

# these graphs represent the change of network between AL and other diets for each tissue
deltagraphs <- list(
  CCR.Blood = d.fgl$CCR$Blood$d.fgl$Delta.graph.connected, 
  ICR.Blood = d.fgl$ICR$Blood$d.fgl$Delta.graph.connected,
  CCR.MFP = d.fgl$CCR$MFP$d.fgl$Delta.graph.connected,
  ICR.MFP = d.fgl$ICR$MFP$d.fgl$Delta.graph.connected
)

targets <- targets.unique # rename
# make new column with ID in KEGG format
targets$gene <- paste("ENTREZID", targets$entrez, sep = ":")
# get all KEGG pathways in iGraph format
pwkegg <- pblapply(pathways("mmusculus", "kegg"),
                   function(x) igraph.from.graphNEL(pathwayGraph(x)))

mir2pw <- function(mir_graph, targets, pw, directed = TRUE) {
  # make sure that graph is directed
  mir_graph <- as.directed(mir_graph, mode = "mutual")
  # keep only target pointed by these miRNA
  targets_of_mirna <- subset(targets, mirna %in% V(mir_graph)$name)

  # build links between miRNA and genes in the pathway
  df <- targets_of_mirna[targets_of_mirna$gene %in% V(pw)$name,
                         c("mirna", "gene")]
  if(nrow(df) == 0) {
    return(list(net = NULL, score = 0))
  }
  mirna_targets_graph <- graph_from_data_frame(df)

  mir_graph |>
    graph.union(mirna_targets_graph) |>
    graph.union(pw) -> g2
  mirnm <- grep("mmu", V(g2)$name)
  genenm <- grep("ENTREZID:", V(g2)$name)

  V(g2)[mirnm]$kind <- "miRNA"
  V(g2)[genenm]$kind <- "gene"

  centrality <- eigen_centrality(g2, directed)
  # centrality is a eigen_centrality object
  score <- centrality$vector[df$gene]
  return(list(net = g2, score = score))
}

all_pathways <- lapply(delta.graphs, function(mir_graph) pblapply(pwkegg,
                         function(pw) mir2pw(mir_graph, targets, pw)))

centralitykegg <- lapply(all_pathways, function(mir_graph) sapply(mir_graph,
                      function(pw) sum(pw$score)))

# df <- pblapply(centralitykegg,
#                function(x)
#                  as.data.frame(x, nm = "score") |>
#                  rownames_to_column(var = "term") |>
#                  arrange(desc(score)))
# 
# mir2kegglist <- list(
#   hif1 = pblapply(deltagraphs, mir2pw, targets = targets, 
#                   pw = pwkegg$`HIF-1 signaling pathway`),
#   mapk = pblapply(deltagraphs, mir2pw, targets = targets, 
#                   pw = pwkegg$`MAPK signaling pathway`),
#   phos = pblapply(deltagraphs, mir2pw, targets = targets, 
#                   pw = pwkegg$`Phosphatidylinositol signaling system`),
#   foxo = pblapply(deltagraphs, mir2pw, targets = targets, 
#                   pw = pwkegg$`FoxO signaling pathway`)
# )
# 
# save(mir2kegglist, file = "Data/mir2kegglist.RData")
# openxlsx::write.xlsx(
#   x = pblapply(centralitykegg, function(x)
#     as.data.frame(x, nm = "score") |> 
#       rownames_to_column(var = "term") |> 
#       arrange(desc(score))),
#   file = "Results/mir2pathway.xlsx"
# )

df0 <- data.frame(CCR.Blood = centralitykegg$CCR.Blood,
                  ICR.Blood = centralitykegg$ICR.Blood,
                  CCR.MFP = centralitykegg$CCR.MFP,
                  ICR.MFP = centralitykegg$ICR.MFP)
openxlsx::write.xlsx(df0, "Results/influence_score.xlsx")

ccr.blood.mfp <- ggplot(data = df0, aes(x = CCR.Blood, y = CCR.MFP)) + 
  geom_jitter() + 
  geom_point(data = df0 |> 
               select(CCR.Blood, CCR.MFP) |> 
               dplyr::filter(CCR.Blood == 0, CCR.MFP > 0),
             aes(x = CCR.Blood, y = CCR.MFP), color = "red") +
  geom_point(data = df0 |> 
               select(CCR.Blood, CCR.MFP) |> 
               dplyr::filter(CCR.Blood > 0, CCR.MFP == 0),
             aes(x = CCR.Blood, y = CCR.MFP), color = "red") + 
  geom_text_repel(aes(label = rownames(df0)), max.overlaps = 15) + 
  # geom_text(aes(label = rownames(df0)), check_overlap = TRUE,) +
  # coord_cartesian(xlim = c(-2, 16)) + 
  labs(title = "Blood vs. MFP (CCR)", 
       subtitle = "Distrubution of Total Influence") +
  theme_bw()

icr.blood.mfp <- ggplot(data = df0, aes(x = ICR.Blood, y = ICR.MFP)) + 
  geom_jitter() +
  geom_point(data = df0 |> 
               select(ICR.Blood, ICR.MFP) |> 
               dplyr::filter(ICR.Blood == 0, ICR.MFP >0),
             aes(x = ICR.Blood, y = ICR.MFP), color = "red") + 
  geom_point(data = df0 |> 
               select(ICR.Blood, ICR.MFP) |> 
               dplyr::filter(ICR.Blood > 0, ICR.MFP == 0),
             aes(x = ICR.Blood, y = ICR.MFP), color = "red") + 
  geom_text_repel(aes(label = rownames(df0)), max.overlaps = 15) +
  labs(title = "Blood vs. MFP (ICR)", 
       subtitle = "Distrubution of Total Influence") +
  theme_bw()

ggsave(filename = "Results/influence_ccr_blood_mfp.png", plot = ccr.blood.mfp, width = 10, height = 10)
ggsave(filename = "Results/influence_icr_blood_mfp.png", plot = icr.blood.mfp, width = 10, height = 10)

ccr.icr.blood <- ggplot(data = df0, aes(x = CCR.Blood, y = ICR.Blood)) + 
  geom_jitter() + 
  geom_point(data = df0 |> 
               select(CCR.Blood, ICR.Blood) |> 
               dplyr::filter(CCR.Blood == 0, ICR.Blood > 0),
             aes(x = CCR.Blood, y = ICR.Blood), color = "red") +
  geom_point(data = df0 |> 
               select(CCR.Blood, ICR.Blood) |> 
               dplyr::filter(CCR.Blood > 0, ICR.Blood == 0),
             aes(x = CCR.Blood, y = ICR.Blood), color = "red") + 
  geom_text_repel(aes(label = rownames(df0)), max.overlaps = 15) + 
  # geom_text(aes(label = rownames(df0)), check_overlap = TRUE,) +
  # coord_cartesian(xlim = c(-2, 16)) + 
  labs(title = "CCR vs ICR (Blood)", 
       subtitle = "Distrubution of Total Influence") +
  theme_bw()

ccr.icr.mfp <- ggplot(data = df0, aes(x = CCR.MFP, y = ICR.MFP)) + 
  geom_jitter() + 
  geom_point(data = df0 |> 
               select(CCR.MFP, ICR.MFP) |> 
               dplyr::filter(CCR.MFP == 0, ICR.MFP > 0),
             aes(x = CCR.MFP, y = ICR.MFP), color = "red") +
  geom_point(data = df0 |> 
               select(CCR.MFP, ICR.MFP) |> 
               dplyr::filter(CCR.MFP > 0, ICR.MFP == 0),
             aes(x = CCR.MFP, y = ICR.MFP), color = "red") + 
  geom_text_repel(aes(label = rownames(df0)), max.overlaps = 15) + 
  # geom_text(aes(label = rownames(df0)), check_overlap = TRUE,) +
  # coord_cartesian(xlim = c(-2, 16)) + 
  labs(title = "CCR vs ICR (MFP)", 
       subtitle = "Distrubution of Total Influence") +
  theme_bw()

ggsave(filename = "Results/influence_ccr_icr_blood.png", plot = ccr.icr.blood, width = 10, height = 10)
ggsave(filename = "Results/influence_ccr_icr_mfp.png", plot = ccr.icr.mfp, width = 10, height = 10)

