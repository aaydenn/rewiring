library(dplyr)
library(biomaRt)

load("Data/d.fgl.RData")
load("Data/targets.Rdata")

delta.centrality <- lapply(list(
  CCR.Blood = d.fgl$CCR$Blood$d.centrality,
  ICR.Blood = d.fgl$ICR$Blood$d.centrality,
  CCR.MFP = d.fgl$CCR$MFP$d.centrality,
  ICR.MFP = d.fgl$ICR$MFP$d.centrality
), function(x)
  x |> mutate(hub = ifelse(knn < degree, 1, 0)) |> dplyr::filter(hub == 1))

eigen.targets <-
  lapply(delta.centrality, function(x)
    targets.unique |> dplyr::filter(mirna %in% (x[x$hub == 1,]$miRNA)))

# openxlsx::write.xlsx(eigen.targets, "Results/eigen.targets.xlsx")

eigen.count <- lapply(eigen.targets, function(x) x |> count(symbol) |> arrange(desc(n)) |> dplyr::filter(n>1))

go <-
  read_delim(
    file = "Data/Databases/gene_association.mgi",
    delim = "\t",
    comment = "!",
    na = "",
    col_names = c(
      "Database Designation",
      "MGI Marker Accession ID",
      "Mouse Marker Symbol",
      "NOT Designation",
      "GO Term ID",
      "MGI Reference Accession ID",
      "GO Evidence Code",
      "Inferred From",
      "Ontology",
      "Mouse Marker Name",
      "Mouse Marker Synonyms (if any)",
      "Mouse Marker Type (gene, transcript, protein)",
      "Taxon",
      "Modification Date",
      "Assigned By",
      "Annotation Extension",	"Gene Product"))

goterm <-
  read_delim(
    file = "Data/Databases/go_terms.mgi",
    delim = "\t",
    col_names = c("Ontology", "GO Term ID", "GO Term")
  )

eigen.go <-
  lapply(eigen.count, function(x)
    go |> dplyr::filter(`Mouse Marker Symbol` %in% x$symbol) |> dplyr::select(`Mouse Marker Symbol`, `GO Term ID`))

