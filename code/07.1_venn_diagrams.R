# make miRNA names from each DE table ####
# "gnames" contains IDs of miRNAs
# load("Data/contrasts.RData")
load("Data/fit.eb.RData")
# names of all de mirnas in a list with the names of contrast
sapply(colnames(fit.eb$contrasts),
       function(x)
         rownames(topTable(
           fit.eb,
           coef = x,
           number = Inf,
           p.value = 0.05
         ))) -> gnames

# diet comparison (fixed age and tissue) ####
# list for venn diagrams ####
venn.diet.genes <- list(
  d.49.Blood = gnames[c("d.AL.49.50.Blood", "d.CCR.49.50.Blood", "d.ICR.49.50.Blood")],
  d.81.Blood = gnames[c("d.AL.81.82.Blood", "d.CCR.81.82.Blood", "d.ICR.81.82.Blood")],
  d.49.Brain = gnames[c("d.AL.49.50.Brain", "d.CCR.49.50.Brain", "d.ICR.49.50.Brain")],
  d.81.Brain = gnames[c("d.AL.81.82.Brain", "d.CCR.81.82.Brain", "d.ICR.81.82.Brain")],
  d.49.MFP   = gnames[c("d.AL.49.50.MFP",   "d.CCR.49.50.MFP",   "d.ICR.49.50.MFP")],
  d.81.MFP   = gnames[c("d.AL.81.82.MFP",   "d.CCR.81.82.MFP",   "d.ICR.81.82.MFP")]
)
# Blood in 49/50
venn.diagram(
  x              = venn.diet.genes$d.49.Blood,
  category.names = c("AL", "CCR", "ICR"),
  main           = "Blood in 49/50",
  filename       = "Results/d.49.Blood.png",
  output         = TRUE,
  euler.d        = FALSE,
  scaled         = FALSE,
  imagetype      = "png" ,
  height         = 1725,
  width          = 1725,
  resolution     = 300,
  cex            = 1.2,
  col            = c("#141542", '#a7a523', '#e6e6c1'),
  fill           = c("#1415427f", '#8a89287f', '#e1e1657f'),
)
# Blood in 81/82
venn.diagram(
  x              = venn.diet.genes$d.81.Blood,
  category.names = c("AL", "CCR", "ICR"),
  main           = "Blood in 81/82",
  filename       = "Results/d.81.Blood.png",
  output         = TRUE,
  euler.d        = FALSE,
  scaled         = FALSE,
  imagetype      = "png" ,
  height         = 1725,
  width          = 1725,
  resolution     = 300,
  cex            = 1.2,
  col            = c("#141542", '#a7a523', '#e6e6c1'),
  fill           = c("#1415427f", '#8a89287f', '#e1e1657f'),
)
# Brain in 49/50
venn.diagram(
  x              = venn.diet.genes$d.49.Brain,
  category.names = c("AL", "CCR", "ICR"),
  main           = "Brain in 49/50",
  filename       = "Results/d.49.Brain.png",
  output         = TRUE,
  euler.d        = FALSE,
  scaled         = FALSE,
  imagetype      = "png" ,
  height         = 1725,
  width          = 1725,
  resolution     = 300,
  cex            = 1.2,
  col            = c("#141542", '#a7a523', '#e6e6c1'),
  fill           = c("#1415427f", '#8a89287f', '#e1e1657f'),
)
# Brain in 81/82
venn.diagram(
  x              = venn.diet.genes$d.81.Brain,
  category.names = c("AL", "CCR", "ICR"),
  main           = "Brain in 81/82",
  filename       = "Results/d.81.Brain.png",
  output         = TRUE,
  euler.d        = FALSE,
  scaled         = FALSE,
  imagetype      = "png" ,
  height         = 1725,
  width          = 1725,
  resolution     = 300,
  cex            = 1.2,
  col            = c("#141542", '#a7a523', '#e6e6c1'),
  fill           = c("#1415427f", '#8a89287f', '#e1e1657f'),
)
# MFP in 49/50
venn.diagram(
  x              = venn.diet.genes$d.49.MFP,
  category.names = c("AL", "CCR", "ICR"),
  main           = "MFP in 49/50",
  filename       = "Results/d.49.MFP.png",
  output         = TRUE,
  euler.d        = FALSE,
  scaled         = FALSE,
  imagetype      = "png" ,
  height         = 1725,
  width          = 1725,
  resolution     = 300,
  cex            = 1.2,
  col            = c("#141542", '#a7a523', '#e6e6c1'),
  fill           = c("#1415427f", '#8a89287f', '#e1e1657f'),
)
# MFP in 81/82
venn.diagram(
  x              = venn.diet.genes$d.81.MFP,
  category.names = c("AL", "CCR", "ICR"),
  main           = "MFP in 81/82",
  filename       = "Results/d.81.MFP.png",
  output         = TRUE,
  euler.d        = FALSE,
  scaled         = FALSE,
  imagetype      = "png" ,
  height         = 1725,
  width          = 1725,
  resolution     = 300,
  cex            = 1.2,
  col            = c("#141542", '#a7a523', '#e6e6c1'),
  fill           = c("#1415427f", '#8a89287f', '#e1e1657f'),
)
# tissue comparison (fixed age and diet) ####
# list for venn diagrams ####
venn.tissue.genes <- list(
  d.AL.49  = gnames[c("d.AL.49.50.Blood" , "d.AL.49.50.Brain" , "d.AL.49.50.MFP")],
  d.AL.81  = gnames[c("d.AL.81.82.Blood" , "d.AL.81.82.Brain" , "d.AL.81.82.MFP")],
  d.CCR.49 = gnames[c("d.CCR.49.50.Blood", "d.CCR.49.50.Brain", "d.CCR.49.50.MFP")],
  d.CCR.81 = gnames[c("d.CCR.81.82.Blood", "d.CCR.81.82.Brain", "d.CCR.81.82.MFP")],
  d.ICR.49 = gnames[c("d.ICR.49.50.Blood", "d.ICR.49.50.Brain", "d.ICR.49.50.MFP")],
  d.ICR.81 = gnames[c("d.ICR.81.82.Blood", "d.ICR.81.82.Brain", "d.ICR.81.82.MFP")]
)
# AL in 49/50
venn.diagram(
  x              = venn.tissue.genes$d.AL.49,
  category.names = c("Blood", "Brain", "MFP"),
  main           = "AL in 49/50",
  filename       = "Results/d.AL.49.png",
  output         = TRUE,
  euler.d        = FALSE,
  scaled         = FALSE,
  imagetype      = "png" ,
  height         = 1725,
  width          = 1725,
  resolution     = 300,
  cex            = 1.2,
  col            = c("#141542", '#a7a523', '#e6e6c1'),
  fill           = c("#1415427f", '#8a89287f', '#e1e1657f'),
)
# AL in 81/82
venn.diagram(
  x              = venn.tissue.genes$d.AL.81,
  category.names = c("Blood", "Brain", "MFP"),
  main           = "AL in 81/82",
  filename       = "Results/d.AL.81.png",
  output         = TRUE,
  euler.d        = FALSE,
  scaled         = FALSE,
  imagetype      = "png" ,
  height         = 1725,
  width          = 1725,
  resolution     = 300,
  cex            = 1.2,
  col            = c("#141542", '#a7a523', '#e6e6c1'),
  fill           = c("#1415427f", '#8a89287f', '#e1e1657f'),
)
# CCR in 49/50
venn.diagram(
  x              = venn.tissue.genes$d.CCR.49,
  category.names = c("Blood", "Brain", "MFP"),
  main           = "CCR in 49/50",
  filename       = "Results/d.CCR.49.png",
  output         = TRUE,
  euler.d        = FALSE,
  scaled         = FALSE,
  imagetype      = "png" ,
  height         = 1725,
  width          = 1725,
  resolution     = 300,
  cex            = 1.2,
  col            = c("#141542", '#a7a523', '#e6e6c1'),
  fill           = c("#1415427f", '#8a89287f', '#e1e1657f'),
)
# CCR in 81/82
venn.diagram(
  x              = venn.tissue.genes$d.CCR.81,
  category.names = c("Blood", "Brain", "MFP"),
  main           = "CCR in 81/82",
  filename       = "Results/d.CCR.81.png",
  output         = TRUE,
  euler.d        = FALSE,
  scaled         = FALSE,
  imagetype      = "png" ,
  height         = 1725,
  width          = 1725,
  resolution     = 300,
  cex            = 1.2,
  col            = c("#141542", '#a7a523', '#e6e6c1'),
  fill           = c("#1415427f", '#8a89287f', '#e1e1657f'),
)
# ICR in 49/50
venn.diagram(
  x              = venn.tissue.genes$d.ICR.49,
  category.names = c("Blood", "Brain", "MFP"),
  main           = "ICR in 49/50",
  filename       = "Results/d.ICR.49.png",
  output         = TRUE,
  euler.d        = FALSE,
  scaled         = FALSE,
  imagetype      = "png" ,
  height         = 1725,
  width          = 1725,
  resolution     = 300,
  cex            = 1.2,
  col            = c("#141542", '#a7a523', '#e6e6c1'),
  fill           = c("#1415427f", '#8a89287f', '#e1e1657f'),
)
# ICR in 81/82
venn.diagram(
  x              = venn.tissue.genes$d.ICR.81,
  category.names = c("Blood", "Brain", "MFP"),
  main           = "ICR in 81/82",
  filename       = "Results/d.ICR.81.png",
  output         = TRUE,
  euler.d        = FALSE,
  scaled         = FALSE,
  imagetype      = "png" ,
  height         = 1725,
  width          = 1725,
  resolution     = 300,
  cex            = 1.2,
  col            = c("#141542", '#a7a523', '#e6e6c1'),
  fill           = c("#1415427f", '#8a89287f', '#e1e1657f'),
)
