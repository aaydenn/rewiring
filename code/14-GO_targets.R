library(dplyr)
library(circlize)


load("result/fgl.RData")
load("result/diff.rank.eigen.RData")
load("result/targets.Rdata")
load("data/go_kegg_list.RData")

#' -----------------------------------
#' @title get KEGG terms of hub miRNAs
#' @description hub mirna enrichment of kegg terms. 
#' miRNAs with degree higher than knn
#' is considered hub


# filter targets for 98 mirna
targets <- targets.unique[targets.unique$mirna %in% colnames(theta[[1]]),]


# list of every genes that targeted by 98 DE mirna
symbol2mir_id <- split(targets$mirna,targets$symbol)


# list of every gene in every KEGG pathway
symbol_in_any_kegg <- KEGG.list.filtered |> unlist() |> unique()


# number of mirnas that targets any gene within any KEGG pathway
num_mirna_total <- symbol2mir_id[symbol_in_any_kegg] |> unlist() |> unique() |> length()


# hypergeometric test
eval_pathway <- function(pw_genes, DE_mirna, symbol2mir_id, num_mirna_total) {
  pw_mirna <- unique(unlist(symbol2mir_id[pw_genes]))
  k <- length(pw_mirna)
  q <- length(intersect(pw_mirna, DE_mirna))
  m <- length(unique(DE_mirna))
  n <- num_mirna_total - m
  p <- phyper(q, m, n, k, lower.tail = FALSE)
  return(c(p.val=p, pw.size=k, de.in.pw=q))
}


# multiple test for every pathway
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


# KEGG enrichment for hub miRNAs in theta graphs
theta.kegg <- lapply(theta.hub, function(x) {pwenrich(x$miRNA) |> filter(p.val != 0)})

openxlsx::write.xlsx(theta.kegg, file = "tables/theta.kegg.xlsx")


# KEGG enrichment for hub miRNAs in delta graphs
delta.kegg <- lapply(delta.hub, function(x) pwenrich(x$miRNA) |> filter(p.val != 0))

openxlsx::write.xlsx(delta.kegg, file = "tables/delta.kegg.xlsx")




# KEGG enrichment for diff.eigen.rank miRNAs in theta graphs

diff.rank.eigen.kegg <- lapply(diff.rank.eigen, function(x) pwenrich(x$miRNA) |> filter(p.val != 0))

openxlsx::write.xlsx(delta.kegg, diff.rank.eigen.kegg,file = "tables/diff.rank.eigen.kegg.xlsx")



save(theta.kegg,delta.kegg,diff.rank.eigen.kegg,file = "result/kegg.enrichment.RData")

#' -----------------------------------
#' @title get GO terms of hub miRNAs
#' @description hub mirna enrichment of go terms. 
#' miRNAs with degree higher than knn
#' is considered hub


eval_go <- function(go_genes, DE_mirna, symbol2mir_id, num_mirna_total) {
  go_mirna <- unique(unlist(symbol2mir_id[go_genes]))
  k <- length(go_mirna)
  q <- length(intersect(go_mirna, DE_mirna))
  m <- length(unique(DE_mirna))
  n <- num_mirna_total - m
  p <- phyper(q, m, n, k, lower.tail = FALSE)
  return(c(p.val=p, go.size=k, de.in.go=q))
}

goenrich <- function(DE_mirna) {
  GO.list.BP |> 
    sapply(eval_go, DE_mirna, symbol2mir_id, num_mirna_total) |> 
    t() |> 
    as.data.frame() -> ans
  ans$GO_Term <- names(GO.list.BP)
  ans$adj.p.val <- p.adjust(ans$p.val, method = "fdr")
  ans <- ans[order(ans$adj.p.val),]
  return(dplyr::select(ans, GO_Term, go.size, de.in.go, p.val, adj.p.val))
}



# GO enrichment for hub miRNAs in theta graphs
theta.go <- lapply(theta.hub, function(x) {goenrich(x$miRNA) |> filter(p.val != 0)})

openxlsx::write.xlsx(theta.go, file = "tables/theta.go.xlsx")


# GO enrichment for hub miRNAs in delta graphs
delta.go <- lapply(delta.hub, function(x) {goenrich(x$miRNA) |> filter(p.val != 0)})

openxlsx::write.xlsx(delta.go, file = "tables/delta.go.xlsx")


# GO enrichment for diff.eigen.rank miRNAs in theta graphs
diff.rank.eigen.go <- lapply(diff.rank.eigen, function(x) goenrich(x$miRNA) |> filter(p.val != 0))

openxlsx::write.xlsx(diff.rank.eigen.go, file = "tables/diff.rank.eigen.go.xlsx")



save(theta.go, delta.go, diff.rank.eigen.go, file = "result/go.enrichment.RData")
