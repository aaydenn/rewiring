library(igraph)
load("Data/d.fgl.Rdata") # first level mirna-mirna (GLASSO) ####
load("Data/targets.Rdata") # second level mirna-gene (mirnatap) ####
# choose only connected graphs
deltagraphs <- list(CCR.Blood = d.fgl$CCR$Blood$d.fgl$Delta.graph.connected, 
                    ICR.Blood = d.fgl$ICR$Blood$d.fgl$Delta.graph.connected,
                    CCR.MFP = d.fgl$CCR$MFP$d.fgl$Delta.graph.connected,
                    ICR.MFP = d.fgl$ICR$MFP$d.fgl$Delta.graph.connected)
# rename
targets.unique |> mutate(mirna = str_remove(mirna, "mmu-")) -> targets

# make graphs targeted mutually
mir_graph <- lapply(deltagraphs, as.directed, mode = "mutual")

# keep only target pointed by these miRNA
mirna_targets_graph <-
  lapply(lapply(mir_graph, function(x)
    subset(targets, mirna %in% V(x)$name)), graph_from_data_frame)

# degree distribution for each DCE
degree_dist <- function(net1, net2, netid, title = names(net1)) {
  # count the frequencies of each degree
  # convert the first column to numbers
  deg.hist1 <-
    as.data.frame(table(degree(net1))) |>
    mutate(Var1 = as.numeric(Var1), id = netid[1])
  
  deg.hist2 <-
    as.data.frame(table(degree(net2))) |>
    mutate(Var1 = as.numeric(Var1), id = netid[2])
  
  deg.hist <- rbind(deg.hist1, deg.hist2)
  
  #plot the degree distribution
  ggplot(deg.hist, aes(x = Var1, y = Freq)) +
    geom_point(aes(color = id)) +
    geom_line(aes(color = id)) +
    scale_x_continuous("Degree k",
                       trans = "log10") +
    scale_y_continuous("Frequency P(k)",
                       trans = "log10") +
    scale_color_manual("Type", values = c("orange", "purple")) +
    labs(title = title, subtitle = "Degree Distribution (log-log)") + 
    theme_bw() + 
    theme(legend.position = c(0.7, 0.9))
}

ggsave(
  filename = "Results/degdist_d.CCR.Blood.png",
  degree_dist(mirna_targets_graph$CCR.Blood, 
              deltagraphs$CCR.Blood, 
              c("miRNA--target", "miRNA--miRNA"), 
              "CCR vs. AL (Blood)"),
  width = 7, height = 7, units = "in", dpi = 320
)
ggsave(
  filename = "Results/degdist_d.ICR.Blood.png",
  degree_dist(mirna_targets_graph$ICR.Blood, 
              deltagraphs$ICR.Blood, 
              c("miRNA--target", "miRNA--miRNA"), 
              "ICR vs. AL (Blood)"),
  width = 7, height = 7, units = "in", dpi = 320
)
ggsave(
  filename = "Results/degdist_d.CCR.MFP.png",
  degree_dist(mirna_targets_graph$CCR.MFP, 
              deltagraphs$CCR.MFP, 
              c("miRNA--target", "miRNA--miRNA"), 
              "CCR vs. AL (MFP)"),
  width = 7, height = 7, units = "in", dpi = 320
)
ggsave(
  filename = "Results/degdist_d.ICR.MFP.png",
  degree_dist(mirna_targets_graph$ICR.MFP, 
              deltagraphs$ICR.MFP, 
              c("miRNA--target", "miRNA--miRNA"), 
              "ICR vs. AL (MFP)"),
  width = 7, height = 7, units = "in", dpi = 320
)

path_dist <- function(net, title) {
  dist <- as.data.frame(table(distances(net)))
  dist[, 1] <- as.numeric(dist[, 1])
  # omit zero length
  dist <- dist[-c(1),]
  # omit infinite length
  dist <- dist[-c(nrow(dist)),]
  ggplot(data = dist, aes(x = Var1, y = log(Freq))) + 
    geom_bar(stat = "identity") +
    labs(title = title, subtitle = "Shortest Path Length Distribution (log)") +
    xlab("Length of shortest path") + 
    ylab("log(Frequency)") +
    theme_bw()
}

ggsave(
  filename = "Results/degdist_d.CCR.Blood.png",
  path_dist(mirna_targets_graph$CCR.Blood, "CCR vs. AL (Blood)"),
  width = 6
)
ggsave(
  filename = "Results/degdist_d.ICR.Blood.png",
  path_dist(mirna_targets_graph$ICR.Blood, "ICR vs. AL (Blood)"),
  width = 6
)
ggsave(
  filename = "Results/degdist_d.CCR.MFP.png",
  path_dist(mirna_targets_graph$CCR.MFP, "CCR vs. AL (MFP)"),
  width = 6
)
ggsave(
  filename = "Results/degdist_d.ICR.MFP.png",
  path_dist(mirna_targets_graph$ICR.MFP, "ICR vs. AL (MFP)"),
  width = 6
)

mir2gene <- function(net, title) {
  sub <- targets |> subset(mirna %in% V(net)$name)
  df <- as.data.frame(table(sub$mirna))
  
  ggplot(df, aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity") + 
    scale_color_manual("Type", values = "brown") +
    labs(title = title) + 
    theme_bw() + 
    theme(legend.position = c(0.9, 0.9))
}

mir2gene(mirna_targets_graph$CCR.Blood, "ccr")

degree_dist_hist <- function(net1, net2, netid, title = names(net1)) {
  # count the frequencies of each degree
  # convert the first column to numbers
  deg.hist1 <-
    as.data.frame(table(degree(net1))) |>
    mutate(Var1 = as.numeric(Var1), id = netid[1])
  
  deg.hist2 <-
    as.data.frame(table(degree(net2))) |>
    mutate(Var1 = as.numeric(Var1), id = netid[2])
  
  deg.hist <- rbind(deg.hist1, deg.hist2)
  
  #plot the degree distribution
  ggplot(deg.hist, aes(x = Var1, y = log(Freq))) +
    geom_col(aes(color = id, fill = id), 
             position = position_dodge(0.8), 
             width = 0.7) +
    labs(title = title, 
         subtitle = "Degree Distribution (log-log)",
         fill = "Interaction") + 
    theme_bw()
}

degree_dist_hist(mirna_targets_graph$CCR.Blood, 
            deltagraphs$CCR.Blood, 
            c("miRNA--target", "miRNA--miRNA"), 
            "CCR vs. AL (Blood)")

degree_dist_hist(mirna_targets_graph$ICR.Blood, 
                 deltagraphs$ICR.Blood, 
                 c("miRNA--target", "miRNA--miRNA"), 
                 "ICR vs. AL (Blood)")

#####

deg.hist1 <-
  

deg.hist2 <-
  

deg.hist <- rbind(as.data.frame(table(degree(deltagraphs$CCR.Blood))) |>
                    mutate(Var1 = as.numeric(Var1), id = "miRNA--miRNA"), 
                  as.data.frame(table(degree(mirna_targets_graph$CCR.Blood))) |>
                    mutate(Var1 = as.numeric(Var1), id = "miRNA--target"))

ggplot(deg.hist, aes(x = Var1, y = Freq)) +
  geom_point(aes(color = id)) +
  geom_line(aes(color = id)) +
  scale_x_continuous("Degree k",
                     trans = "log10") +
  scale_y_continuous("Frequency P(k)",
                     trans = "log10") +
  scale_color_manual("Type", values = c("orange", "purple")) +
  labs(title = title, subtitle = "Degree Distribution (log-log)") + 
  theme_bw() + 
  theme(legend.position = c(0.7, 0.9))
