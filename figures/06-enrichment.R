library(ggplot2)
library(dplyr)
library(stringr)

load("result/fgl.RData")
load("result/diff.rank.eigen.RData")
load("result/targets.Rdata")
load("data/go_kegg_list.RData")
load("result/kegg.enrichment.RData")
load("result/go.enrichment.RData")

#' Theta GO plots

theta.go.top <- lapply(theta.go, function(x) 
  x |> arrange(adj.p.val) |> head(10) |> 
    setNames(c("GO_Term", "GO_Size","miRNA_Size","Pval","Adj_Pval"))) 


theta.go.plot <- lapply(theta.go.top, function(x) {
  ggplot(x,aes(reorder(GO_Term,-log10(Adj_Pval)),-log10(Adj_Pval),fill=miRNA_Size))+
    geom_bar(stat="identity")+
    theme_bw()+
    coord_flip()+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 20))+
    labs(x="",y="")
  })

p <- (theta.go.plot$Blood.AL + labs(title = "Blood", subtitle = "AL", x="GO_Term")  | theta.go.plot$MFP.AL + labs(title = "MFP", subtitle = "AL")) / 
  (theta.go.plot$Blood.CCR + labs(subtitle = "CCR", x="GO_Term") | theta.go.plot$MFP.CCR + labs(subtitle = "CCR")) / 
  (theta.go.plot$Blood.ICR + labs(subtitle = "ICR", y="log(Adj_Pval)", x="GO_Term") | theta.go.plot$MFP.ICR + labs(subtitle = "ICR", y="log(Adj_Pval)"))


ggsave(filename = "figures/theta.go.png", plot = p, width = 8, height = 15, units = "in",dpi = 300)



#' Theta KEGG plots

theta.kegg.top <- lapply(theta.kegg, function(x) 
  x |> arrange(adj.p.val) |> head(10) |> 
    setNames(c("KEGG_Term", "KEGG_Size","miRNA_Size","Pval","Adj_Pval"))) 


theta.kegg.plot <- lapply(theta.kegg.top, function(x) {
  ggplot(x,aes(reorder(KEGG_Term,-log10(Adj_Pval)),-log10(Adj_Pval),fill=miRNA_Size))+
    geom_bar(stat="identity")+
    scale_fill_gradient(low = "salmon1", high = "salmon4")+
    theme_bw()+
    coord_flip()+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 20))+
    labs(x="", y="")
})

p1 <- (theta.kegg.plot$Blood.AL + labs(title = "Blood", subtitle = "AL", x="KEGG_Term")  | theta.kegg.plot$MFP.AL + labs(title = "MFP", subtitle = "AL")) / 
  (theta.kegg.plot$Blood.CCR + labs(subtitle = "CCR", x="KEGG_Term") | theta.kegg.plot$MFP.CCR + labs(subtitle = "CCR")) / 
  (theta.kegg.plot$Blood.ICR + labs(subtitle = "ICR", y="log(Adj_Pval)", x="KEGG_Term") | theta.kegg.plot$MFP.ICR + labs(subtitle = "ICR", y="log(Adj_Pval)"))


ggsave(filename = "figures/theta.kegg.png", plot = p1, width = 8, height = 15, units = "in",dpi = 300)



#' Diff eigen.rank GO plots
 
diff.rank.eigen.go.top <- lapply(diff.rank.eigen.go, function(x) 
  x |> arrange(adj.p.val) |> head(10) |> 
    setNames(c("GO_Term", "GO_Size","miRNA_Size","Pval","Adj_Pval"))) 


diff.rank.eigen.go.plot <- lapply(diff.rank.eigen.go.top, function(x) {
  ggplot(x,aes(reorder(GO_Term,-log10(Adj_Pval)),-log10(Adj_Pval),fill=miRNA_Size))+
    geom_bar(stat="identity")+
    scale_fill_gradient(low = "slateblue1", high = "slateblue4")+
    theme_bw()+
    coord_flip()+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 20))+
    labs(x="",y="")
})

p2 <- (diff.rank.eigen.go.plot$diff.rank.Blood.CCR.AL + labs(title = "Blood", subtitle = "CCR AL", x="GO_Term") | 
         diff.rank.eigen.go.plot$diff.rank.MFP.CCR.AL + labs(title = "MFP", subtitle = "CCR AL")) / 
  (diff.rank.eigen.go.plot$diff.rank.Blood.ICR.AL + labs(subtitle = "ICR AL", y="log(Adj_Pval)", x="GO_Term") | 
     diff.rank.eigen.go.plot$diff.rank.MFP.ICR.AL + labs(subtitle = "ICR AL", y="log(Adj_Pval)"))


ggsave(filename = "figures/diff.rank.go.png", plot = p2, width = 8, height = 15, units = "in",dpi = 300)


#' Diff. eigen.rank KEGG plots


diff.rank.eigen.kegg.top <- lapply(diff.rank.eigen.kegg, function(x) 
  x |> arrange(adj.p.val) |> head(10) |> 
    setNames(c("KEGG_Term", "KEGG_Size","miRNA_Size","Pval","Adj_Pval"))) 


diff.rank.eigen.kegg.plot <- lapply(diff.rank.eigen.kegg.top, function(x) {
  ggplot(x,aes(reorder(KEGG_Term,-log10(Adj_Pval)),-log10(Adj_Pval),fill=miRNA_Size))+
    geom_bar(stat="identity")+
    scale_fill_gradient(low = "sienna1", high = "sienna4")+
    theme_bw()+
    coord_flip()+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 20))+
    labs(x="", y="")
})

p3 <- (diff.rank.eigen.kegg.plot$diff.rank.Blood.CCR.AL + labs(title = "Blood", subtitle = "CCR AL", x="KEGG_Term") | 
         diff.rank.eigen.kegg.plot$diff.rank.MFP.CCR.AL + labs(title = "MFP", subtitle = "CCR AL")) / 
  (diff.rank.eigen.kegg.plot$diff.rank.Blood.ICR.AL + labs(subtitle = "ICR AL", y="log(Adj_Pval)", x="KEGG_Term") | 
     diff.rank.eigen.kegg.plot$diff.rank.MFP.ICR.AL + labs(subtitle = "ICR AL", y="log(Adj_Pval)"))

ggsave(filename = "figures/diff.rank.kegg.png", plot = p3, width = 8, height = 15, units = "in",dpi = 300)

