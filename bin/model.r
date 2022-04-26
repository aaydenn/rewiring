#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library("optparse")))
suppressWarnings(suppressMessages(library("limma")))
suppressWarnings(suppressMessages(library("dplyr")))
suppressWarnings(suppressMessages(library("readr")))

option_list <- list(
    make_option(
        c("-n", "--normalized-data"), type = "character", 
        default = NULL, metavar = "/path/to/rdata/",
        help = "R data matrix object that contains normalized microarray data"),
    make_option(
        c("-f", "--features"), type = "character", 
        default = NULL, metavar = "/path/to/rdata/",
        help = "R data matrix object that contains features of the samples"),
    make_option(
        c("-x", "--contrasts"), type="character", 
        default = NULL, metavar = "/path/to/contrast_matrix/",
        help = "Contrast matrix for experimental design")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

load(opt$n)
load(opt$f)

design <- model.matrix(
    ~ 0 + diet:age:tissue, 
    data = features
)

design <- design[, colSums(design)>0]

colnames(design) |> 
  str_remove("diet") |>
  str_remove("tissue") |>
  str_remove("age") |> 
  str_replace_all("/","_") |> 
  str_replace_all(":","_") -> colnames(design)

cat("model matrix created!\n")

fit <- lmFit(mouse_data, design)

if(!is.null(opt$x)) {
    contrasts <- makeContrasts(
        levels = design,
        as.matrix(opt$contrasts))
} else {
    contrasts <- makeContrasts(
        levels = design,
        d_AL_49_50_Blood  = AL_49_50_Blood  - BASELINE_10_Blood,
        d_CCR_49_50_Blood = CCR_49_50_Blood - BASELINE_10_Blood,
        d_ICR_49_50_Blood = ICR_49_50_Blood - BASELINE_10_Blood,
        d_AL_81_82_Blood  = AL_81_82_Blood  - BASELINE_10_Blood,
        d_CCR_81_82_Blood = CCR_81_82_Blood - BASELINE_10_Blood,
        d_ICR_81_82_Blood = ICR_81_82_Blood - BASELINE_10_Blood,
        d_AL_49_50_MFP    = AL_49_50_MFP    - BASELINE_10_MFP,
        d_CCR_49_50_MFP   = CCR_49_50_MFP   - BASELINE_10_MFP,
        d_ICR_49_50_MFP   = ICR_49_50_MFP   - BASELINE_10_MFP,
        d_AL_81_82_MFP    = AL_81_82_MFP    - BASELINE_10_MFP,
        d_CCR_81_82_MFP   = CCR_81_82_MFP   - BASELINE_10_MFP,
        d_ICR_81_82_MFP   = ICR_81_82_MFP   - BASELINE_10_MFP
    )
}


fit_cont <- contrasts.fit(fit, contrasts)
fit_eb <- eBayes(fit_cont)

cat("fitted to the contrast!\n")

ft_eb_table <- sapply(
    colnames(fit_eb$contrasts), function(i) 
        rownames_to_column(topTable(fit_eb, coef = i,
            number = Inf, p.value = 0.05), "miRNA"))

for (i in 1:length(fit_eb_table)) {
    if (nrow(fit_eb_table[[i]])>0) {
        write_tsv(fit_eb_table[[i]], names(fit_eb_table[i]), escape = "none")
    } else {
        cat("No significant miRNA found: ", names(fit_eb_table[i]))
    }
}

sig_ids <- sapply(
    colnames(fit_eb$contrasts), function(x) 
        rownames(topTable(fit_eb, coef = x, 
            number = Inf, p.value = 0.05))
)

sig_ids |> 
    unlist() |> 
    unique() -> de_list

save(de_list, file = "de_list.rda")

cat("significant miRNAs found!\n")