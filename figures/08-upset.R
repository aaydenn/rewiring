# make miRNA names from each DE table ####
# "gnames" contains IDs of miRNAs
# load("Data/contrasts.RData")

library(UpSetR)
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

blood <- venn_table("d.AL.49.50.Blood", "d.AL.81.82.Blood",
                    "d.CCR.49.50.Blood", "d.CCR.81.82.Blood",
                    "d.ICR.49.50.Blood", "d.ICR.81.82.Blood")

colnames(blood) <- c("miRNA",
                     " AL.50", " AL.81",
                     "CCR.50", "CCR.81",
                     "ICR.50", "ICR.81")

mfp <- venn_table("d.AL.49.50.MFP", "d.AL.81.82.MFP",
                  "d.CCR.49.50.MFP", "d.CCR.81.82.MFP",
                  "d.ICR.49.50.MFP", "d.ICR.81.82.MFP")

colnames(mfp) <- c("miRNA",
                   " AL.50", " AL.81",
                   "CCR.50", "CCR.81",
                   "ICR.50", "ICR.81")

png(filename = "figures/blood_upset.png", width = 8, height = 5, units = "in", res = 300)
upset(blood, 
      sets = c(" AL.50", "CCR.50", "ICR.50",
               " AL.81", "CCR.81", "ICR.81"), 
      keep.order = TRUE, 
      mainbar.y.label = "Intersection (Blood)", 
      sets.x.label = NULL,
      matrix.color = "gray11", point.size = 5, text.scale = 1.35)
dev.off()

png(filename = "figures/mfp_upset.png", width = 8, height = 5, units = "in", res = 300)
upset(mfp, 
      sets = c(" AL.50", "CCR.50", "ICR.50",
               " AL.81", "CCR.81", "ICR.81"), 
      keep.order = TRUE, 
      mainbar.y.label = "Intersection (MFP)", 
      sets.x.label = NULL,
      matrix.color = "gray11", point.size = 5, text.scale = 1.35)
dev.off()

openxlsx::write.xlsx(blood,file = "tables/blood.upset.matrix.xlsx")
openxlsx::write.xlsx(mfp,file = "tables/mfp.upset.matrix.xlsx")
