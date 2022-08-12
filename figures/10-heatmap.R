library(pheatmap)
library(viridis)


load("result/mouse-data-features.RData")
load("result/DE_miRNAs.RData")
# colnames(mouse.Data) |> 
#   stringr::str_remove(".CEL") |> 
#   stringr::str_remove(".cel") |> 
#   stringr::str_remove(".ga") -> colnames(mouse.Data)


colnames(mouse.Data) <- toupper(colnames(mouse.Data))


rownames(features) <- colnames(mouse.Data)
# rownames(features)[rownames(features)=="T71Blood"] <- "ICRR_T71_Blood_w80"
# rownames(features)[rownames(features)=="T25W80Blood"] <- "ICRRF_T25_Blood_w80"


# modify diet column of features
df <- features |> 
        mutate(diet = ifelse(diet == "ICRR" | diet == "ICRRF", "ICR", diet),
               age = ifelse(age == "49/50", "50", ifelse(age == "81/82", "81", "10"))) |> 
        select(diet, tissue, age) |> 
        filter(tissue!="Brain") |> 
  rename(Tissue = tissue, Diet = diet, Age = age)


mouse.Data |>
  as.data.frame() |> 
  filter(rownames(mouse.Data) %in% c(DE.blood, DE.mfp)) -> arrays

arrays <- arrays[features$tissue != "Brain"]
rownames(arrays) <- stringr::str_remove(rownames(arrays), "mmu-")

colour = list(
  Tissue = c(Blood = "brown4", MFP = "darkolivegreen"),
  Age = c(`10` = "honeydew", `50` = "cadetblue3", `81` = "darkslategray"),
  Diet = c(BASELINE = "plum1", AL = "plum3", CCR = "violetred1", ICR = "violetred4")
)

png(filename = "figures/heatmap0.png", 4, 8, units = "in", res = 320)
pheatmap(arrays, show_colnames = F,
         show_rownames = F,
         annotation_col = df, 
         annotation_colors = colour, border_color = NA)
dev.off()






# kegg heatmap

load("result/theta.influence.RData")

pw <- do.call(cbind, theta.influence) |> as.data.frame() |> filter(Blood.AL > quantile(Blood.AL, 0.90))



png(filename = "figures/pw.png", 8, 8, units = "in", res = 320)

pheatmap(pw, angle_col = 45, border_color = NA)

dev.off()
