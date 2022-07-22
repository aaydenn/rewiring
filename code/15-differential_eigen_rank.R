load("result/fgl.RData")
load("result/go.enrichment.RData")
load("result/kegg.enrichment.RData")
load("result/theta.influence.RData")

#' ------------------------------------------------------------------
#' condition dependent differential eigen centralitty ranks of miRNAs

theta.centrality.ranked <- lapply(theta.centrality, function(x) x |> mutate(eigen.rank=1:nrow(x)) |> select(miRNA,eigen.rank))

# Blood

diff.rank.theta.go.Blood <- inner_join(theta.go.ranked$Blood.AL,theta.centrality.ranked$Blood.CCR,by="miRNA") |> 
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

openxlsx::write.xlsx(diff.rank.eigen, "tables/diff.rank.eigen.xlsx")



#' -------------------------------------------------------
#' condition dependent differential enriched GO term ranks

theta.go.ranked <- lapply(theta.go, function(x) x |> mutate(adj.p.val.rank=1:nrow(x)))



#' Blood

diff.rank.theta.go.Blood.CCR.AL <- inner_join(theta.go.ranked$Blood.AL,theta.go.ranked$Blood.CCR,by="GO_Term") |> 
  mutate(d.rank = adj.p.val.rank.x - adj.p.val.rank.y) |> 
  mutate(z.rank = (d.rank- mean(d.rank))/sd(d.rank)) |> filter(abs(z.rank)>1)

diff.rank.theta.go.Blood.ICR.AL <- inner_join(theta.go.ranked$Blood.AL,theta.go.ranked$Blood.ICR,by="GO_Term") |> 
  mutate(d.rank = adj.p.val.rank.x - adj.p.val.rank.y) |> 
  mutate(z.rank = (d.rank- mean(d.rank))/sd(d.rank)) |> filter(abs(z.rank)>1)


#' MFP

diff.rank.theta.go.MFP.CCR.AL <- inner_join(theta.go.ranked$MFP.AL,theta.go.ranked$MFP.CCR,by="GO_Term") |> 
  mutate(d.rank = adj.p.val.rank.x - adj.p.val.rank.y) |> 
  mutate(z.rank = (d.rank- mean(d.rank))/sd(d.rank)) |> filter(abs(z.rank)>1)

diff.rank.theta.go.MFP.ICR.AL <- inner_join(theta.go.ranked$MFP.AL,theta.go.ranked$MFP.ICR,by="GO_Term") |> 
  mutate(d.rank = adj.p.val.rank.x - adj.p.val.rank.y) |> 
  mutate(z.rank = (d.rank- mean(d.rank))/sd(d.rank)) |> filter(abs(z.rank)>1)

diff.rank.theta.go <- list(
  theta.go.Blood.CCR.AL = diff.rank.theta.go.Blood.CCR.AL,
  theta.go.Blood.ICR.AL = diff.rank.theta.go.Blood.ICR.AL,
  theta.go.MFP.CCR.AL = diff.rank.theta.go.MFP.CCR.AL,
  theta.go.MFP.ICR.AL = diff.rank.theta.go.MFP.ICR.AL
)

save(diff.rank.theta.go, file = "result/diff.rank.theta.go.RData")

openxlsx::write.xlsx(diff.rank.theta.go, "tables/diff.rank.theta.go.xlsx")



#' -------------------------------------------------------
#' condition dependent differential enriched KEGG term ranks

theta.kegg.ranked <- lapply(theta.kegg, function(x) x |> mutate(adj.p.val.rank=1:nrow(x)) |> select(term,adj.p.val.rank))



#' Blood

diff.rank.theta.kegg.Blood.CCR.AL <- inner_join(theta.kegg.ranked$Blood.AL,theta.kegg.ranked$Blood.CCR,by="term") |> 
  mutate(d.rank = adj.p.val.rank.x - adj.p.val.rank.y) |> 
  mutate(z.rank = (d.rank- mean(d.rank))/sd(d.rank)) |> filter(abs(z.rank)>1)

diff.rank.theta.kegg.Blood.ICR.AL <- inner_join(theta.kegg.ranked$Blood.AL,theta.kegg.ranked$Blood.ICR,by="term") |> 
  mutate(d.rank = adj.p.val.rank.x - adj.p.val.rank.y) |> 
  mutate(z.rank = (d.rank- mean(d.rank))/sd(d.rank)) |> filter(abs(z.rank)>1)


#' MFP

diff.rank.theta.kegg.MFP.CCR.AL <- inner_join(theta.kegg.ranked$MFP.AL,theta.kegg.ranked$MFP.CCR,by="term") |> 
  mutate(d.rank = adj.p.val.rank.x - adj.p.val.rank.y) |> 
  mutate(z.rank = (d.rank- mean(d.rank))/sd(d.rank)) |> filter(abs(z.rank)>1)

diff.rank.theta.kegg.MFP.ICR.AL <- inner_join(theta.kegg.ranked$MFP.AL,theta.kegg.ranked$MFP.ICR,by="term") |> 
  mutate(d.rank = adj.p.val.rank.x - adj.p.val.rank.y) |> 
  mutate(z.rank = (d.rank- mean(d.rank))/sd(d.rank)) |> filter(abs(z.rank)>1)

diff.rank.theta.kegg <- list(
  theta.kegg.Blood.CCR.AL = diff.rank.theta.kegg.Blood.CCR.AL,
  theta.kegg.Blood.ICR.AL = diff.rank.theta.kegg.Blood.ICR.AL,
  theta.kegg.MFP.CCR.AL = diff.rank.theta.kegg.MFP.CCR.AL,
  theta.kegg.MFP.ICR.AL = diff.rank.theta.kegg.MFP.ICR.AL
)

save(diff.rank.theta.kegg, file = "result/diff.rank.theta.kegg.RData")

openxlsx::write.xlsx(diff.rank.theta.kegg, "tables/diff.rank.theta.kegg.xlsx")


#' ------------------------------------------------
#' condition dependent "influence" on KEGG term ranks
 

theta.influence <- pblapply(theta.influence, function(x)
  as.data.frame(x, nm = "score") |>
    rownames_to_column(var = "term") |>
    arrange(desc(score)))


theta.influence.ranked <- lapply(theta.influence, function(x) x |> mutate(score.rank=1:nrow(x)) |> select(term, score.rank))

#' Blood

diff.rank.influence.Blood.CCR.AL <- inner_join(theta.influence.ranked$Blood.AL,theta.influence.ranked$Blood.CCR,by="term") |> 
  mutate(d.rank = score.rank.x - score.rank.y) |> 
  mutate(z.rank = (d.rank- mean(d.rank))/sd(d.rank)) |> filter(abs(z.rank)>1)

diff.rank.influence.Blood.ICR.AL <- inner_join(theta.influence.ranked$Blood.AL,theta.influence.ranked$Blood.ICR,by="term") |> 
  mutate(d.rank = score.rank.x - score.rank.y) |> 
  mutate(z.rank = (d.rank- mean(d.rank))/sd(d.rank)) |> filter(abs(z.rank)>1)


#' MFP

diff.rank.influence.MFP.CCR.AL <- inner_join(theta.influence.ranked$MFP.AL,theta.influence.ranked$MFP.CCR,by="term") |> 
  mutate(d.rank = score.rank.x - score.rank.y) |> 
  mutate(z.rank = (d.rank- mean(d.rank))/sd(d.rank)) |> filter(abs(z.rank)>1)

diff.rank.influence.MFP.ICR.AL <- inner_join(theta.influence.ranked$MFP.AL,theta.influence.ranked$MFP.ICR,by="term") |> 
  mutate(d.rank = score.rank.x - score.rank.y) |> 
  mutate(z.rank = (d.rank- mean(d.rank))/sd(d.rank)) |> filter(abs(z.rank)>1)

diff.rank.theta.influence <- list(
  influence.Blood.CCR.AL = diff.rank.influence.Blood.CCR.AL,
  influence.Blood.ICR.AL = diff.rank.influence.Blood.ICR.AL,
  influence.MFP.CCR.AL = diff.rank.influence.MFP.CCR.AL,
  influence.MFP.ICR.AL = diff.rank.influence.MFP.ICR.AL
)

save(diff.rank.theta.influence, file = "result/diff.rank.theta.influence.RData")

openxlsx::write.xlsx(diff.rank.theta.influence, "tables/diff.rank.theta.influence.xlsx")
