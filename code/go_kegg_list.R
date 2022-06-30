#' -----------------------------
#' @title read GO terms from MGI

go <- read_delim(file = "data/gene_association.mgi.gz", 
                 delim = "\t", comment = "!", na = "",
                 col_names = c("Database_Designation","MGI_Marker_Accession_ID","Mouse_Marker_Symbol","NOT Designation","GO_Term_ID","MGI_Reference_Accession_ID","GO_Evidence_Code","Inferred_From","Ontology","Mouse_Marker_Name","Mouse_Marker_Synonyms","Mouse_Marker_Type","Taxon","Modification_Date","Assigned_By","Annotation_Extension","Gene_Product")) |> 
  mutate(Ontology = ifelse(Ontology == "P", "BP",
                           ifelse(Ontology == "F", "MF", "CC")))

goterm <- read_delim(file = "data/go_terms.mgi", 
                     delim = "\t", col_names = c("Ontology", "GO_Term_ID", "GO_Term")) |> 
  mutate(Ontology = ifelse(Ontology == "Biological Process", "BP",
                           ifelse(Ontology == "Molecular Function", "MF", "CC")))

GO.table <- go |> 
  left_join(goterm, by = c("GO_Term_ID", "Ontology"))

GO.table.BP <- GO.table |> filter(Ontology == "BP")


GO.list <- lapply(split(x=GO.table$Mouse_Marker_Symbol, f=GO.table$GO_Term), unique)

GO.list.BP <- lapply(split(x=GO.table.BP$Mouse_Marker_Symbol, f=GO.table.BP$GO_Term), unique)


rm(go,goterm)

#' ----------------------------------------------
#' @title read KEGG terms from 'graphite' package

library(graphite)

pw <- convertIdentifiers(pathways("mmusculus", "kegg"), "symbol")

pwlist <- lapply(pw, function(x) nodes(x) |> str_remove("SYMBOL:"))

KEGG.list <- pwlist

discard <- c(grep("metabolism",names(KEGG.list),value = T),
             grep("disease",names(KEGG.list),value = T), 
             grep("infection",names(KEGG.list),value = T),
             names(which(lengths(KEGG.list)==0)))

KEGG.list.filtered <- KEGG.list[names(KEGG.list) %in% discard == FALSE]

save(GO.table,GO.list,GO.table.BP,GO.list.BP,
     KEGG.list,KEGG.list.filtered,file = "data/go_kegg_list.RData")
