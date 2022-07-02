library(dplyr)
library(circlize)

load("result/jgl.RData")
load("result/diff.rank.eigen.RData")
load("result/targets.Rdata")
load("data/go_kegg_list.RData")

#' -----------------------------------
#' @title get KEGG terms of hub miRNAs
#' @description hub mirna enrichment of kegg terms. 
#' miRNAs with degree higher than knn
#' is considered hub


symbol2mir_id <- split(targets.unique$mirna,targets.unique$symbol)
symbol_in_any_kegg <- KEGG.list.filtered |> unlist() |> unique()
num_mirna_total <- symbol2mir_id[symbol_in_any_kegg] |> unlist() |> unique() |> length()

eval_pathway <- function(pw_genes, DE_mirna, symbol2mir_id, num_mirna_total) {
  pw_mirna <- unique(unlist(symbol2mir_id[pw_genes]))
  k <- length(pw_mirna)
  q <- length(intersect(pw_mirna, DE_mirna))
  m <- length(unique(DE_mirna))
  n <- num_mirna_total - m
  p <- phyper(q, m, n, k, lower.tail = FALSE)
  return(c(p.val=p, pw.size=k, de.in.pw=q))
}

pwenrich <- function(DE_mirna) {
  KEGG.list.filtered |> 
    sapply(eval_pathway, DE_mirna, symbol2mir_id, num_mirna_total) |> 
    t() |> 
    as.data.frame() -> ans
  ans$term <- names(KEGG.list.filtered)
  ans$adj.p.val <- p.adjust(ans$p.val, method = "fdr")
  ans <- ans[order(ans$adj.p.val),]
  return(dplyr::select(ans, term, pw.size, de.in.pw, p.val, adj.p.val))
}


theta.hub <- lapply(theta.centrality, function(x) x |> filter(hub == 1))

delta.hub <- lapply(delta.centrality, function(x) x |> filter(hub == 1))



theta.kegg <- lapply(theta.hub, function(x) {pwenrich(x$miRNA) |> filter(p.val != 0)})

openxlsx::write.xlsx(theta.kegg, file = "tables/theta.kegg.xlsx")



delta.kegg <- lapply(delta.hub, function(x) pwenrich(x$miRNA) |> filter(p.val != 0))

openxlsx::write.xlsx(delta.kegg, file = "tables/delta.kegg.xlsx")



save(theta.kegg,delta.kegg,file = "result/kegg.enrichment.RData")

#' -----------------------------------
#' @title get KEGG terms of hub miRNAs
#' @description hub mirna enrichment of go terms. 
#' miRNAs with degree higher than knn
#' is considered hub


eval_go <- function(go_genes, DE_mirna, symbol2mir_id, num_mirna_total) {
  pw_mirna <- unique(unlist(symbol2mir_id[go_genes]))
  k <- length(pw_mirna)
  q <- length(intersect(pw_mirna, DE_mirna))
  m <- length(unique(DE_mirna))
  n <- num_mirna_total - m
  p <- phyper(q, m, n, k, lower.tail = FALSE)
  return(c(p.val=p, pw.size=k, de.in.pw=q))
}

goenrich <- function(DE_mirna) {
  GO.list.BP |> 
    sapply(eval_go, DE_mirna, symbol2mir_id, num_mirna_total) |> 
    t() |> 
    as.data.frame() -> ans
  ans$GO_Term <- names(GO.list.BP)
  ans$adj.p.val <- p.adjust(ans$p.val, method = "fdr")
  ans <- ans[order(ans$adj.p.val),]
  return(dplyr::select(ans, GO_Term, pw.size, de.in.pw, p.val, adj.p.val))
}




theta.go <- lapply(theta.hub, function(x) {goenrich(x$miRNA) |> filter(p.val != 0)})

openxlsx::write.xlsx(theta.go, file = "tables/theta.go.xlsx")




delta.go <- lapply(delta.hub, function(x) {goenrich(x$miRNA) |> filter(p.val != 0)})

openxlsx::write.xlsx(delta.go, file = "tables/delta.go.xlsx")




save(theta.go, delta.go, file = "result/go.enrichment.RData")

















# filter genes targeted more than 1
# 90% = 2  95% = 2 100% = 5

hub.targets <- data.frame(symbol=c(Arid1="Arid1",Hif1a="Hif1a",Hif1a="Hif1a",Tet2="Tet2"))

hub.count <- lapply(hub.targets, function(x)
  x |> 
    count(symbol) |> 
    arrange(desc(n)) |> 
    dplyr::filter(n > 2))


# get annotations of frequent target genes

hub.go <- lapply(hub.count, function(x)
  GO |> 
    dplyr::filter(Mouse_Marker_Symbol %in% x$symbol) |> 
    dplyr::select(Mouse_Marker_Symbol, GO_Term_ID, GO_Term, Ontology) |> 
    distinct() |> 
    filter(Ontology == "BP"))

go.count <- 
  lapply(hub.go, function(x)
    x |> count(GO_Term_ID) |>
      filter(n > 1))

go.table <-
  lapply(setNames(names(hub.go), names(hub.go)), function(x)
    hub.go[[x]] |>
      filter(GO_Term_ID %in% go.count[[x]]$GO_Term_ID) |>
      with(table(GO_Term_ID, Mouse_Marker_Symbol)))



