library(org.Mm.eg.db)
library(miRNAtap)
library(annotate)
library(pbapply)

load("Data/DE_miRNAs.RData")

# TODO: describe logic behind gene selection

## `targets` is a data.frame/matrix
## We ignore the last 2 columns (why?)
geom_mean_of_sample <- function(targets) {
  n <- ncol(targets) - 2 # numb. of sources - matrix w/ score of different sources and rank (last two columns are rank_product and rank final)
  if (n>1) { # does not work for targets with only one source
    sampling <- sapply(1:n, function(j) sample(targets[,j], # for every source, simulate random target list
                                           size = 1000, replace = TRUE))
    weight <- rowSums((!is.na(sampling))) # sum of non-NA rows
    valid_row <- weight > 0 # exclude rows with zero weight
    weight <- weight[valid_row]
    sampling <- sampling[valid_row, ]
    geom_meon <- exp(rowMeans(log(sampling), na.rm = TRUE))
    return(geom_meon/weight) #geometric means corrected to weight
  } else
    return(NULL) # simulate sampling for a single column
}

permutat_threshold <- function(targets, signif = 0.05) {
  gms <- geom_mean_of_sample(targets = targets)
  thr <- sort(gms)[length(gms) * signif] # threshold value for sampled targets
  valid <- targets[, ncol(targets) - 1] <= thr
  return(targets[valid,]) # return targets above threshold
}

get_targets_list <- function(miRNAs) {
  ## in 'getPredictedTargets' function, a line checks whether it is
  ## a mirna or not if it has less then three letters. I propose this:
  ## *if(nchar(mirma)<3|is.na(mirna))*
  ## in same function it takes mirna string between 5 and 15
  ## which is not the case for all mirnas. I propose this:
  ## *mirna = gsub('-3p|-5p', '', mirna)*
  # targets genes for each mirna. A list of vectors/matrices
  list_of_targets <- pblapply(miRNAs, function(miRNA) {
    data.frame(getPredictedTargets(miRNA, species = "mmu"))
  })
  names(list_of_targets) <- miRNAs
  list_of_targets <- pblapply(list_of_targets, permutat_threshold) # implicit 5% significance
  for (miRNA in miRNAs) {
    if (nrow(list_of_targets[[miRNA]])>0) {
      list_of_targets[[miRNA]] <- data.frame(
        mirna = miRNA, entrez = rownames(list_of_targets[[miRNA]]),
        symbol = annotate::getSYMBOL(rownames(list_of_targets[[miRNA]]), 'org.Mm.eg.db')
      )
    } else {
      warning("No target for: ", names(list_of_targets[miRNA]))
    }
  }
  return(list_of_targets)
}

# make one long data frame with all target genes for all miRNA
predicted_targets_df <- function(miRNAs) {
  df <- data.frame(do.call(rbind, get_targets_list(miRNAs)))
  if(nrow(df)>0) {
    rownames(df) <- NULL
    df$type <- "predicted"
  }
  return(df)
}

# this section finds validated targets of miRNAs in a vector
validated_df <- function(miRNAs) {
  mti <- readRDS("Data/mirtarbase_mmu.rds") # database is miRTarBase
  cols <- c("miRNA", "Target Gene (Entrez ID)", "Target Gene")
  df <- unique(mti[mti$miRNA %in% miRNAs, cols])
  if(nrow(df)>0) {
    colnames(df) <- c("mirna","entrez","symbol")
    df$type <- "validated"
  }
  return(df)
}

#####
predicted_targets.unique <- predicted_targets_df(gnames.unique)
validated_targets.unique <- validated_df(gnames.unique)
targets.unique <- rbind(predicted_targets.unique, validated_targets.unique)

#####
predicted_targets_list  <- lapply(gnames, predicted_targets_df)
validated_targets_list <- lapply(gnames, validated_df)
targets_list <- mapply(rbind, predicted_targets_list, validated_targets_list, SIMPLIFY = FALSE)

save(targets.unique, targets_list, file = "Data/targets.Rdata")

#####
# get_targets_list <- function(miRNAs) {
#   require(org.Mm.eg.db)
#   require(crayon)
#   require(annotate)
#   require(miRNAtap)
#   tryCatch({
#     if (length(miRNAs) == 0) {
#       return(NULL)
#     }
# 
#     genes <- list()
#     for (i in 1:length(miRNAs)) {
#       targets <- getPredictedTargets(mirna = miRNAs[i], species ="mmu",
#                                      min_src = 2)
#       targets <- permutat_threshold(targets) # implicit 5% significance
#       if (length(targets) == 0) {
#         cat(red("No targets for ", miRNAs[i], "\n"))
#       } else {
#         genes[[i]] <- targets
#         names(genes)[[i]] <- miRNAs[i]
#         cat(length(genes[[i]]), "target genes for miRNA: ", miRNAs[i], "\n")
#       }
#     }
#     #to do: do these with lapply!
#     table <- list()
#     for (j in 1:length(genes)) {
#       table[[j]] <- cbind(rep(names(genes)[j], length(rownames(genes[[j]]))),
#                           rownames(genes[[j]]))
#     } 
#     result <- do.call("rbind.data.frame", table)
#     names(result) <- c("miRNAId", "EntrezId")
#     symbol <- getSYMBOL(result$EntrezId, 'org.Mm.eg.db')
#     if (length(symbol) == 0) {
#       cat(red("No mappings for target genes"))
#       return(result)
#     } else {
#       result$Symbol <- symbol
#     }
#     cat(blue("Mapping genes finished\n"))
#     return(result)
#   }, error=function(e){
#     cat(red("ERROR in provided list:", conditionMessage(e), "\n"))})
# }
#####