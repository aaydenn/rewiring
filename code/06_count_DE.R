load("Data/fit.eb.RData")
## function for a nice presentation of DE miRNAs
## adjusted p value is 0.05
# de <- function(fit, ind, p = 0.05) {d <- list()
#   for (i in seq_along(ind)) {d[[i]] <- topTable(fit, coef = unlist(ind[i]),number = Inf,p.value = p)
#   }
#   return(setNames(d, ind))
# }
# DE.table <- # store DE mirnas in a list of data frames
#   setNames(de(fit.eb, ind = colnames(fit.eb$contrasts)), colnames(fit.eb$contrasts))
DE.table <- sapply(colnames(fit.eb$contrasts), function(i)
    rownames_to_column(topTable(fit.eb, coef = i,
                       number = Inf, p.value = 0.05), "miRNA"))

# DE.table <- lapply(DE.table, function(x) rownames_to_column(x, "miRNA"))
# DE.table.up <- lapply(DE.table, function(x) x |> filter(x$logFC > 0)) # upregulated mirnas
# DE.table.dn <- lapply(DE.table, function(x) x |> filter(x$logFC < 0)) # downregulated mirnas
data.frame(
  num.ALL = sapply(DE.table, nrow), # we count the total number of DE miRNA for p.value 0.05
  num.DN = sapply(DE.table, function(x) x |> filter(x$logFC < 0) |> nrow()), # for every contrast
                  # we count the number of DE miRNA for p.value 0.05
                  # and fold change smaller than 0
  num.UP = sapply(DE.table, function(x) x |> filter(x$logFC > 0) |> nrow()), # for every contrast
                  # we count the number of DE miRNA for p.value 0.05
                  # and fold change larger than 0
  num.DN.lfc1 = sapply(DE.table, function(x) x |> filter(x$logFC < -1) |> nrow()), # for every contrast
                  # we count the number of DE miRNA for p.value 0.05
                  # and absolute fold change larger than 1
                  # and filter fold change smaller than 1
  num.UP.lfc1 = sapply(DE.table, function(x) x |> filter(x$logFC > 1) |> nrow())
                  # for every contrast we count the number of DE miRNA for p.value 0.05
                  # absolute fold change larger than 1 and filter fold change larger than 1
) -> DE.number

gnames <- sapply(colnames(fit.eb$contrasts), function(x)
  rownames(topTable(fit.eb, coef = x,number = Inf, p.value = 0.05)))

gnames |> unlist() |> unique() -> gnames.unique

gnames[c(
    "d.AL.49.50.Blood","d.CCR.49.50.Blood",
    "d.ICR.49.50.Blood","d.AL.81.82.Blood",
    "d.CCR.81.82.Blood","d.ICR.81.82.Blood")] |> unlist() |> unique() -> DE.blood

gnames[c(
    "d.AL.49.50.Brain","d.CCR.49.50.Brain",
    "d.ICR.49.50.Brain","d.AL.81.82.Brain",
    "d.CCR.81.82.Brain","d.ICR.81.82.Brain")] |> unlist() |> unique() -> DE.brain

gnames[c(
    "d.AL.49.50.MFP","d.CCR.49.50.MFP",
    "d.ICR.49.50.MFP","d.AL.81.82.MFP",
    "d.CCR.81.82.MFP","d.ICR.81.82.MFP")] |> unlist() |> unique() -> DE.mfp

save(DE.table, DE.number, gnames, gnames.unique, DE.blood, DE.brain, DE.mfp, file = "Data/DE_miRNAs.RData")
