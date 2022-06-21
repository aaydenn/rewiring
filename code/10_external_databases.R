# make mirna-gene-tf graph ####
library(org.Mm.eg.db)
library(miRBaseConverter)

mppi <- read_delim("Data/Databases/mippie_ppi_v1_0.tsv", 
# third level gene-gene. mippie curates mouse protein protein interaction from various sources ####
                   delim = "\t",col_names = T)

mppi |> # select only "BioGRID" and "MINT"
  filter(provider == "BioGrid" | provider == "MINT") |> # select only "BioGRID" and "MINT"
  mutate(source = annotate::getSYMBOL(as.character(entrezA), "org.Mm.eg.db"), 
         target = annotate::getSYMBOL(as.character(entrezB), "org.Mm.eg.db"), 
         interaction = "interacts") |> 
  dplyr::select(source, target, interaction) |> unique() -> mppi.edge

tibble(id = c(mppi.edge$source, mppi.edge$target),
                    group = "gene") |> unique() -> mppi.node

read_delim("Data/Databases/Mus_musculus_TF.txt", delim = "\t") |> 
  # names of every mouse TFs "http://bioinfo.life.hust.edu.cn/AnimalTFDB/"
  dplyr::select(Symbol, `Entrez ID`) -> tfbd 

transmir <- read_delim("Data/Databases/mmu_TransmiR.tsv", delim = "\t", 
# fourth level miRNA-TF. transmir curates tf -> mirna regulation ####
                       col_names = c("TF name", "miRNA name", "TSS",
                                     "Binding site", "Action type", "SRAID/PMID",
                                     "Evidence", "Tissue", "Species"))
transmir |> 
  dplyr::mutate(source = `TF name`, 
                target = miRNA_PrecursorToMature(`miRNA name`)$Mature1,
                interaction = ifelse(`Action type` == "Regulation", "interacts", 
                                     ifelse(`Action type` == "Regulation(feedback)", "interacts", 
                                            ifelse(`Action type` == "Activation", "activates", 
                                                   ifelse(`Action type` == "Activation(feedback)", "activates", 
                                                          ifelse(`Action type` == "Repression", "inhibits", "inhibits")))))) |> 
  dplyr::select(source,target,interaction) |> unique() -> transmir.edge

tibble(id = c(transmir.edge$source, transmir.edge$target), 
       group = c(rep("TF", nrow(transmir.edge)), 
                 rep("miRNA", nrow(transmir.edge)))) |> unique() -> transmir.node

trrust <- read_delim("Data/Databases/trrust_rawdata.mouse.tsv", delim = "\t", 
# fifth level TF-gene. TRRUST curates tf - gene regulation ####
           col_names = c("Interactor 1", "Interactor 2", 
                         "Regulation", "References (PMID)"))

trrust |> dplyr::filter(Regulation != "Unknown") |> 
  dplyr::mutate(source = `Interactor 1`,
                target = `Interactor 2`,
                interaction = ifelse(Regulation == "Activation", "activates", "inhibits")) |>
  dplyr::select(source, target, interaction) |> unique() -> trrust.edge

tibble(id = c(trrust.edge$source, trrust.edge$target), 
       group = c(rep("TF", nrow(trrust.edge)), 
                 rep("gene", nrow(trrust.edge)))) |> unique() -> trrust.node

save(mppi, mppi.edge, mppi.node, transmir, transmir.edge, transmir.node,
     trrust, trrust.edge,  trrust.node, file = "Data/databases.RData")


