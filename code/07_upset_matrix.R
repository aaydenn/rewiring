# make miRNA names from each DE table ####
# "gnames" contains IDs of miRNAs
# load("Data/contrasts.RData")
library(UpSetR)
load("Data/DE_miRNAs.RData")
load("Data/mouse-data-features.RData")
# names of all de mirnas in a list with the names of contrast
sapply(colnames(fit.eb$contrasts),
       function(x)
         rownames(topTable(
           fit.eb,
           coef = x,
           number = Inf,
           p.value = 0.05
         ))) -> gnames

# functions to retrieve venn-like tables ####
# with this table we know which miRNA belongs to which group
## table for two contrasts
venn_table_dou <-
  function(c1, c2, gn = gnames, mouse_annot = mouse_4.1) {
    names <- unique(unlist(gn[c(c1, c2)]))
    transcr <- mouse_annot$`Transcript ID(Array Design)`
    names(transcr) <- mouse_annot$`Transcript ID(Array Design)`
    
    data.frame(
      Transcr = transcr[names],
      C1 = names %in% gn[[c1]],
      C2 = names %in% gn[[c2]],
      row.names = NULL
    )
  }

##table for three contrasts
venn_table_tri <-
  function(c1, c2, c3, gn = gnames, mouse_annot = mouse_4.1) {
    names <- unique(unlist(gn[c(c1, c2, c3)]))
    transcr <- mouse_annot$`Transcript ID(Array Design)`
    names(transcr) <- mouse_annot$`Transcript ID(Array Design)`

    data.frame(
      Transcr = transcr[names],
      C1 = names %in% gn[[c1]],
      C2 = names %in% gn[[c2]],
      C3 = names %in% gn[[c3]],
      row.names = NULL
    )
  }
# tables for desired comparisons ####
# comparisons for diet: age and tissue fixed
venn.diet <- list(
  d.49.Blood = venn_table_tri("d.AL.49.50.Blood", "d.CCR.49.50.Blood", "d.ICR.49.50.Blood"),
  d.81.Blood = venn_table_tri("d.AL.81.82.Blood", "d.CCR.81.82.Blood", "d.ICR.81.82.Blood"),
  d.49.Brain = venn_table_tri("d.AL.49.50.Brain", "d.CCR.49.50.Brain", "d.ICR.49.50.Brain"),
  d.81.Brain = venn_table_tri("d.AL.81.82.Brain", "d.CCR.81.82.Brain", "d.ICR.81.82.Brain"),
  d.49.MFP   = venn_table_tri("d.AL.49.50.MFP"  , "d.CCR.49.50.MFP"  , "d.ICR.49.50.MFP"),
  d.81.MFP   = venn_table_tri("d.AL.81.82.MFP"  , "d.CCR.81.82.MFP"  , "d.ICR.81.82.MFP")
)
# comparisons for tissue: age and diet fixed
venn.tissue <- list(
  d.AL.49  = venn_table_tri("d.AL.49.50.Blood" , "d.AL.49.50.Brain" , "d.AL.49.50.MFP"),
  d.AL.81  = venn_table_tri("d.AL.81.82.Blood" , "d.AL.81.82.Brain" , "d.AL.81.82.MFP"),
  d.CCR.49 = venn_table_tri("d.CCR.49.50.Blood", "d.CCR.49.50.Brain", "d.CCR.49.50.MFP"),
  d.CCR.81 = venn_table_tri("d.CCR.81.82.Blood", "d.CCR.81.82.Brain", "d.CCR.81.82.MFP"),
  d.ICR.49 = venn_table_tri("d.ICR.49.50.Blood", "d.ICR.49.50.Brain", "d.ICR.49.50.MFP"),
  d.ICR.81 = venn_table_tri("d.ICR.81.82.Blood", "d.ICR.81.82.Brain", "d.ICR.81.82.MFP")
)
# save tables to rdata
save(venn.diet, venn.tissue, file = "Results/venn_tables.RData")


##table for all contrasts (upset)
venn_table <-
  function(c1, c2, c3, c4,c5, c6, gn = gnames, mouse_annot = mouse_4.1) {
    names <- unique(unlist(gn[c(c1, c2, c3)]))
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
                     "AL.50",
                     "CCR.50",
                     "ICR.50",
                     "AL.80",
                     "CCR.80",
                     "ICR.80")

mfp <- venn_table("d.AL.49.50.MFP", "d.AL.81.82.MFP",
                  "d.CCR.49.50.MFP", "d.CCR.81.82.MFP",
                  "d.ICR.49.50.MFP", "d.ICR.81.82.MFP")

colnames(mfp) <- c("miRNA",
                   "AL.50",
                   "CCR.50",
                   "ICR.50",
                   "AL.80",
                   "CCR.80",
                   "ICR.80")

png(filename = "blood_upset.png", width = 10, height = 10, units = "in", res = 320)
upset(blood, 
      sets = c("AL.50", "CCR.50", "ICR.50",
               "AL.80", "CCR.80", "ICR.80"), 
      keep.order = TRUE, 
      mainbar.y.label = "Intersection (Blood)", 
      sets.x.label = NULL,
      matrix.color = "gray11", point.size = 5, text.scale = 3)
dev.off()

png(filename = "mfp_upset.png", width = 10, height = 10, units = "in", res = 320)
upset(mfp, 
      sets = c("AL.50", "CCR.50", "ICR.50",
               "AL.80", "CCR.80", "ICR.80"), 
      keep.order = TRUE, 
      mainbar.y.label = "Intersection (MFP)", 
      sets.x.label = NULL,
      matrix.color = "gray11", point.size = 5, text.scale = 3)
dev.off()
