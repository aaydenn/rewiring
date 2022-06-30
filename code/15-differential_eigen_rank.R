tmp <- theta.centrality$Blood.AL

tmp |> mutate(rank=1:nrow(tmp))

theta.centrality.ranked <- lapply(theta.centrality, function(x) x |> mutate(eigen.rank=1:nrow(x)) |> select(miRNA,eigen.rank))

#' ------------------------------------------------------------------
#' condition dependent differential eigen centralitty ranks of miRNAs


#' Blood

diff.rank.Blood.CCR.AL <- inner_join(theta.centrality.ranked$Blood.AL,theta.centrality.ranked$Blood.CCR,by="miRNA") |> 
  mutate(d.rank = eigen.rank.x - eigen.rank.y) |> 
  mutate(z.rank = (d.rank- mean(d.rank))/sd(d.rank)) |> filter(abs(z.rank)>1)

diff.rank.Blood.ICR.AL <- inner_join(theta.centrality.ranked$Blood.AL,theta.centrality.ranked$Blood.ICR,by="miRNA") |> 
  mutate(d.rank = eigen.rank.x - eigen.rank.y) |> 
  mutate(z.rank = (d.rank- mean(d.rank))/sd(d.rank)) |> filter(abs(z.rank)>1)


#' MFP
 
diff.rank.MFP.CCR.AL <- inner_join(theta.centrality.ranked$MFP.AL,theta.centrality.ranked$MFP.CCR,by="miRNA") |> 
  mutate(d.rank = eigen.rank.x - eigen.rank.y) |> 
  mutate(z.rank = (d.rank- mean(d.rank))/sd(d.rank)) |> filter(abs(z.rank)>1)

diff.rank.MFP.ICR.AL <- inner_join(theta.centrality.ranked$MFP.AL,theta.centrality.ranked$MFP.ICR,by="miRNA") |> 
  mutate(d.rank = eigen.rank.x - eigen.rank.y) |> 
  mutate(z.rank = (d.rank- mean(d.rank))/sd(d.rank)) |> filter(abs(z.rank)>1)

diff.rank.eigen <- list(
    diff.rank.Blood.CCR.AL = diff.rank.Blood.CCR.AL,
    diff.rank.Blood.ICR.AL = diff.rank.Blood.ICR.AL,
    diff.rank.MFP.CCR.AL = diff.rank.MFP.CCR.AL,
    diff.rank.MFP.ICR.AL = diff.rank.MFP.ICR.AL
  )

save(diff.rank.eigen, file = "result/diff.rank.eigen.RData")

openxlsx::write.xlsx(diff.rank.eigen, "tables/diff.eigen.rank.xlsx")
