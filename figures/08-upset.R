# make miRNA names from each DE table ####
# "gnames" contains IDs of miRNAs
# load("Data/contrasts.RData")

library(ComplexUpset)
load("result/DE_miRNAs.RData")
load("result/mouse-data-features.RData")


##table for all contrasts (upset)
venn_table <-
  function(c1, c2, c3, c4,c5, c6, gn = gnames, mouse_annot = mouse_4.1) {
    names <- unique(unlist(gn[c(c1, c2, c3, c4,c5, c6)]))
    transcr <- mouse_annot$`Transcript ID(Array Design)`
    names(transcr) <- mouse_annot$`Transcript ID(Array Design)`
    
    data.frame(
      Transcr = transcr[names],
      C1 = as.integer(names %in% gn[[c1]]),
      C2 = as.integer(names %in% gn[[c2]]),
      C3 = as.integer(names %in% gn[[c3]]),
      C4 = as.integer(names %in% gn[[c4]]),
      C5 = as.integer(names %in% gn[[c5]]),
      C6 = as.integer(names %in% gn[[c6]]),
      row.names = NULL
    )
  }



# Blood

blood <- venn_table("d.AL.49.50.Blood", "d.AL.81.82.Blood",
                    "d.CCR.49.50.Blood", "d.CCR.81.82.Blood",
                    "d.ICR.49.50.Blood", "d.ICR.81.82.Blood")

colnames(blood) <- c("miRNA",
                     " AL.50", " AL.81",
                     "CCR.50", "CCR.81",
                     "ICR.50", "ICR.81")


b <- do.call(rbind,DE.table[1:6]) |> 
  mutate(feat = rep(names(DE.table[1:6]),sapply(DE.table[1:6],nrow))) |> 
  mutate(Regulation = ifelse(logFC>0, "up","down")) |> 
  select(miRNA,Regulation) |> 
  inner_join(blood, by="miRNA") |> 
  distinct()



png(filename = "figures/blood_upset.png", width = 8, height = 5, units = "in", res = 300)

upset(data = b,
  intersect = colnames(b[-c(1,2)]),
  name='Dietary Groups',
  base_annotations = list('Number of miRNAs' = intersection_size(counts = FALSE, mapping = aes(fill = Regulation)) + 
                            scale_fill_manual(values=c('up'='salmon2', 'down'='royalblue'))),
  sort_sets = FALSE,
  width_ratio = 0.1
)

dev.off()




# MFP

mfp <- venn_table("d.AL.49.50.MFP", "d.AL.81.82.MFP",
                  "d.CCR.49.50.MFP", "d.CCR.81.82.MFP",
                  "d.ICR.49.50.MFP", "d.ICR.81.82.MFP")

colnames(mfp) <- c("miRNA",
                   " AL.50", " AL.81",
                   "CCR.50", "CCR.81",
                   "ICR.50", "ICR.81")

m <- do.call(rbind, DE.table[13:18]) |> 
  mutate(feat = rep(names(DE.table[13:18]),sapply(DE.table[13:18],nrow))) |> 
  mutate(Regulation = ifelse(logFC > 0, "up","down")) |> 
  select(miRNA, Regulation) |> 
  inner_join(mfp, by="miRNA") |> 
  distinct()



png(filename = "figures/mfp_upset.png", width = 8, height = 5, units = "in", res = 300)

upset(data = m,
      intersect = colnames(m[-c(1,2)]),
      name='Dietary Groups',
      base_annotations = list('Number of miRNAs' = intersection_size(counts = FALSE, mapping = aes(fill = Regulation)) + 
                                scale_fill_manual(values=c('up'='salmon2', 'down'='royalblue'))),
      sort_sets = FALSE,
      width_ratio = 0.1)

dev.off()

openxlsx::write.xlsx(blood,file = "tables/blood.upset.matrix.xlsx")
openxlsx::write.xlsx(mfp,file = "tables/mfp.upset.matrix.xlsx")
