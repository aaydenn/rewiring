# sqlite.path="https://dev.mirnet.ca/miRNet/resources/sqlite/"
# mir2gene <- read_delim("Data/Databases/mmu-mir2gene.tsv", quote = "")
# mir2gene$mir_id <- gsub("r","R",mir2gene$mir_id)
load("Data/fit.eb.RData")
load("Data/targets.Rdata")
load("Data/kegg_mmu.rda")
library(pbapply)
# names of all de mirnas in a list with the names of contrast
sapply(colnames(contrasts), 
       function(x) rownames(topTable(fit.eb, coef = x, 
                                     number = Inf, p.value = 0.05))) -> gnames
# object tlif contains predicted and validated 
# targets for every miRNA
entrez2mir_id <- split(mir2gene$mir_id, mir2gene$entrez)
# any genes in all KEGG pathways
kegg$sets |> unlist() |> unique() -> entrez_in_any_kegg
# 
entrez2mir_id[entrez_in_any_kegg] |>
    unlist() |> unique() |> length() -> num_mirna_total
# pathway enrichment ####
# eval_pathway parameters:
# pw_genes:        entrez ids in given KEGG pathway as list (kegg$sets)
# de_mirna:        differentially expressed mirnas as vector
# entrez2mir_id:   entrez ids in each mirna as list
# num_mirna_total: number of mirnas with the targets genes which belongs any KEGG 
eval_pathway <- function(pw_genes, DE_mirna, entrez2mir_id, num_mirna_total) {
    pw_mirna <- unique(unlist(entrez2mir_id[pw_genes]))
    k <- length(pw_mirna)
    q <- length(intersect(pw_mirna, DE_mirna))
    m <- length(unique(DE_mirna))
    n <- num_mirna_total - m
    p <- phyper(q, m, n, k, lower.tail = FALSE)
    return(c(p.val=p, pw.size=k, de.in.pw=q))
}

pwenrich <- function(DE_mirna) {
    kegg$sets |> 
        sapply(eval_pathway, 
               DE_mirna, 
               entrez2mir_id, 
               num_mirna_total) |> 
        t() |> 
        as.data.frame() -> ans
    ans$term <- kegg$term
    ans$adj.p.val <- p.adjust(ans$p.val, method = "fdr")
    ans <- ans[order(ans$adj.p.val),]
    return(dplyr::select(ans, term, pw.size, de.in.pw, p.val, adj.p.val))
}
kegg.enrichment <- pblapply(gnames, function(x) pwenrich(unlist(x)))

save(kegg.enrichment, file ="kegg.enrichment.RData")

# visualization
enrichviz <- function(df, num = 10, ord = "p.val", title = "df") {
    df[1:as.integer(num),] |> 
        arrange(ord) |> 
        ggplot(aes(x = reorder(`term`, -`p.val`), y = -log10(p.val))) + 
        geom_col() + 
        coord_flip() + 
        ggtitle(title)
}

p <- lapply(kegg.enrichment, function(x) enrichviz(x, title = "enriched pathways"))

ggsave(names(p)[1], p[[1]], "png", "Results", dpi = "retina")