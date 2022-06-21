library(pheatmap)

load("Data/mouse-data-features.RData")
load("Data/DE_miRNAs.RData")
# colnames(mouse.Data) |> 
#   stringr::str_remove(".CEL") |> 
#   stringr::str_remove(".cel") |> 
#   stringr::str_remove(".ga") -> colnames(mouse.Data)

colnames(mouse.Data) <- toupper(colnames(mouse.Data))
rownames(features) <- colnames(mouse.Data)
# rownames(features)[rownames(features)=="T71Blood"] <- "ICRR_T71_Blood_w80"
# rownames(features)[rownames(features)=="T25W80Blood"] <- "ICRRF_T25_Blood_w80"

# modify diet column of features
features$diet <- ifelse(features$diet == "ICRR",  "ICR", features$diet)
features$diet <- ifelse(features$diet == "ICRRF", "ICR", features$diet)

mouse.Data |>
  as.data.frame() |> 
  dplyr::filter(rownames(mouse.Data) %in% c(DE.blood, DE.mfp)) -> arrays

arrays <- arrays[features$tissue != "Brain"]
rownames(arrays) <- stringr::str_remove(rownames(arrays), "mmu-")

colour = list(
  tissue = c(Blood = "#cb410b", MFP = "#82e382"),
  age = c(`10` = "#dd779e", `49/50` = "#f2a682", `81/82` = "#b5ffb2"),
  diet = c(BASELINE = "#dd779e", AL = "#23e016", CCR = "#5d0aaf", ICR = "#4dcfd6")
)

png(filename = "Results/heatmap.png", 10, 8, units = "in", res = 320)
pheatmap(t(arrays), show_rownames = F,
         annotation_row = features[c("diet","tissue")], 
         annotation_colors = colour,
         angle_col = 45,
         fontsize_col = 8)
dev.off()
