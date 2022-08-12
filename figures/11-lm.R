# setwd("/Volumes/GoogleDrive/.shortcut-targets-by-id/1bvKREGKXkHjg7zHGepx0DMnvoFfIoFdf/CR-miRNA/Results/")

library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)


# calculate linear regression for the rank/eigen between diets

get_intreval <- function(data, v1, v2) {
  ans <- data.frame(v = "data", id = data$id, x = data[[v1]], y = data[[v2]])
  ans <- na.omit(ans)
  mdl <- lm(y ~ x, data = ans)
  extra <- 0.05 * diff(range(ans$x))
  newdata <- seq(from=min(ans$x)-extra, to=max(ans$x)+extra, length=20)
  newdata <- data.frame(x=newdata)
  interv <- as.data.frame(predict.lm(mdl, newdata, interval = "predict"))
  interv2 <- as.data.frame(predict.lm(mdl, ans, interval = "predict"))
  # table(interv2$lwr > data$CCR.rank, interv2$upr < data$CCR.rank)
  return(list(ans = cbind(ans, interv2),
              interv = rbind(ans, data.frame(v="lwr", id=NA, x=newdata[[1]], y=interv$lwr),
                           data.frame(v="upr", id=NA, x=newdata[[1]], y=interv$upr))))
}




#' ----------------
#' miRNA regression
#' above or below confidence interval means re-ranked miRNAs



# read centrality table

blood <- read_excel("rankinghubs.xlsx", skip = 2)

mfp <- read_excel("rankinghubs.xlsx", sheet = "MFP", skip = 2)

labels <- c("id", 
  paste(rep(c("AL", "CCR", "ICR"), each=3),
    rep(c("rank", "hub", "eigen"), times=3), 
      sep="."))

colnames(blood) <- colnames(mfp) <- labels



blood$AL.eigen[is.na(blood$AL.hub)] <- NA
blood$CCR.eigen[is.na(blood$CCR.hub)] <- NA
blood$ICR.eigen[is.na(blood$ICR.hub)] <- NA

blood$id |> stringr::str_remove("mmu-") -> blood$id



mfp$AL.eigen[is.na(mfp$AL.hub)] <- NA
mfp$CCR.eigen[is.na(mfp$CCR.hub)] <- NA
mfp$ICR.eigen[is.na(mfp$ICR.hub)] <- NA

mfp$id |> stringr::str_remove("mmu-") -> mfp$id

# plot(CCR.eigen ~ AL.eigen, blood)
# m <- lm(CCR.eigen ~ AL.eigen, blood)
# abline(m)
# newdata <- data.frame(AL.eigen = seq(from = -0.05, to = 1.05, by = 0.05))
# interv <- as.data.frame(predict.lm(m, newdata, interval = "predict"))
# lines(newdata$AL.eigen, interv$lwr, col="red")
# lines(newdata$AL.eigen, interv$upr, col="red")
# interv2 <- as.data.frame(predict.lm(m, interval = "predict"))
# table(interv2$lwr > blood$CCR.eigen, interv2$upr < blood$CCR.eigen)

# base plot

a <- get_intreval(blood, "AL.eigen", "CCR.eigen")

plot(y~x, data=a$ans, subset = v=="data", xlab = "AL.eigen", ylab = "CCR.eigen")
lines(y~x, data=a$interv, subset = v=="lwr", col="red")
lines(y~x, data=a$interv, subset = v=="upr", col="red")

a$ans |> filter(y<lwr | y>upr)



# ggplot
# x = CCR, y = AL

eigeninterval <- list(Blood.CCR.AL  = get_intreval(blood, "CCR.eigen", "AL.eigen")$ans,
                      Blood.ICR.AL  = get_intreval(blood, "ICR.eigen", "AL.eigen")$ans,
                      Blood.ICR.CCR = get_intreval(blood, "ICR.eigen", "CCR.eigen")$ans,
                      MFP.CCR.AL  = get_intreval(mfp, "CCR.eigen", "AL.eigen")$ans,
                      MFP.ICR.AL  = get_intreval(mfp, "ICR.eigen", "AL.eigen")$ans,
                      MFP.ICR.CCR = get_intreval(mfp, "ICR.eigen", "CCR.eigen")$ans)

rankinterval <- list(Blood.CCR.AL  = get_intreval(blood, "CCR.rank", "AL.rank")$ans,
                     Blood.ICR.AL  = get_intreval(blood, "ICR.rank", "AL.rank")$ans,
                     Blood.ICR.CCR = get_intreval(blood, "ICR.rank", "CCR.rank")$ans,
                     MFP.CCR.AL  = get_intreval(mfp, "CCR.rank", "AL.rank")$ans,
                     MFP.ICR.AL  = get_intreval(mfp, "ICR.rank", "AL.rank")$ans,
                     MFP.ICR.CCR = get_intreval(mfp, "ICR.rank", "CCR.rank")$ans)


# compare eigen centrality

eigenpl <- lapply(eigeninterval, function(i) {
  ggplot(i) + 
    geom_smooth(aes(x = x, y = y), method = "lm", se = FALSE, color = "darkgrey") +
    # geom_line(aes(x = x, y = lwr)) + 
    # geom_line(aes(x = x, y = upr)) + 
    geom_ribbon(aes(x = x,y = y, ymin = upr, ymax = lwr), fill = alpha("grey",0.5)) +
    geom_point(aes(x = x, y = y)) + 
    geom_point(aes(x = x, y = y, size = 1.1), 
               filter(i, y<lwr | y>upr), col = alpha("firebrick",0.5)) + 
    geom_text_repel(aes(x = x, y = y, label = id), 
               filter(i, y<lwr | y>upr), col = "firebrick") + 
    theme_bw() + 
    theme(legend.position = "NA", 
          text = element_text(size = 15))
})


p0 <- (((eigenpl$Blood.CCR.AL + labs(x="",y="eigen (AL)")) + 
         (eigenpl$Blood.ICR.AL + labs(x="",y=""))) / 
        ((eigenpl$MFP.CCR.AL + labs(x="eigen (CCR)",y="eigen (AL)")) + 
           (eigenpl$MFP.ICR.AL + labs(x="eigen (ICR)",y="")))) + 
  plot_annotation(tag_levels = "A")

ggsave(p0, filename = "eigen_manuscript.svg", device = "svg", width = 8,height = 8,dpi = 320)




# compare rank

rankpl <- lapply(rankinterval, function(i) {
  ggplot(i) + 
    geom_smooth(aes(x = x, y = y), method = "lm", se = FALSE, color = "darkgrey") +
    # geom_line(aes(x = x, y = lwr)) + 
    # geom_line(aes(x = x, y = upr)) + 
    geom_ribbon(aes(x = x,y = y, ymin = upr, ymax = lwr), fill = alpha("grey",0.5)) +
    geom_point(aes(x = x, y = y)) + 
    geom_point(aes(x = x, y = y, size = 1.1), 
               filter(i, y<lwr | y>upr), col = alpha("firebrick",0.5)) + 
    geom_text_repel(aes(x = x, y = y, label = id), 
                    filter(i, y<lwr | y>upr), col = "firebrick") + 
    
    theme_bw() + 
    theme(legend.position = "NA", 
          text = element_text(size = 15))
})

p <- (((rankpl$Blood.CCR.AL + labs(x="",y="rank (AL)")) + 
    (rankpl$Blood.ICR.AL + labs(x="",y=""))) / 
    ((rankpl$MFP.CCR.AL + labs(x="rank (CCR)",y="rank (AL)")) + 
       (rankpl$MFP.ICR.AL + labs(x="rank (ICR)",y="")))) + 
  plot_annotation(tag_levels = "A")

ggsave(p, filename = "rank_manuscript.svg", device = "svg", width = 8,height = 8,dpi = 320)





#' ----------------
#' influence (pw) regression
#' above or below confidence interval means re-ranked pws


pwys <- read_excel("theta.influence.xlsx", sheet = "All", skip = 1)
colnames(pwys) <- c("id", paste(rep(c("blood", "mfp"), each=4),
                                rep(c("AL", "CCR", "ICR", "avg"), times=2), sep="."))

a <- get_intreval(pwys, "blood.AL", "blood.CCR")
a$ans |> filter(y<lwr | y>upr) |> select(id, x, y)

a <- get_intreval(pwys, "blood.AL", "blood.ICR")
a$ans |> filter(y<lwr | y>upr) |> select(id, x, y)

a <- get_intreval(pwys, "mfp.AL", "mfp.CCR")
a$ans |> filter(y<lwr | y>upr) |> select(id, x, y)

a <- get_intreval(pwys, "mfp.AL", "mfp.ICR")
a$ans |> filter(y<lwr | y>upr) |> select(id, x, y)

