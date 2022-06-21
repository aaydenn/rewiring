library(readr)
library(tidyverse)

mouse <- read_delim("Data/samplelist.tsv")

for (i in 1:ncol(mouse)) {
  mouse[[i]][is.na(mouse[[i]])] <- 0
}

mouse <- subset(mouse, sacrifice_week > 0)

t.test(breast_cancer ~ diet,
       data = mouse,
       subset = (diet == "AL" | diet == "CCR"))
t.test(breast_cancer ~ diet,
       data = mouse,
       subset = (diet == "AL" | diet == "ICR-R"))
t.test(breast_cancer ~ diet,
       data = mouse,
       subset = (diet == "AL" | diet == "ICR-RF"))

mouse$CR <-
  ifelse(mouse$diet == "CCR" |
           mouse$diet == "ICR-R" | mouse$diet == "ICR-RF",
         "any_CR",
         mouse$diet)

t.test(
  breast_cancer ~ CR,
  data = mouse,
  alternative = "greater",
  subset = diet != "Baseline"
)
t.test(
  breast_cancer ~ CR,
  data = mouse,
  alternative = "greater",
  subset = (diet != "Baseline" & sacrifice_week %in% c(49, 50))
)
t.test(
  breast_cancer ~ CR,
  data = mouse,
  alternative = "greater",
  subset = (diet != "Baseline" & sacrifice_week %in% c(81, 82))
)
t.test(breast_cancer ~ CR,
       data = mouse,
       subset = (diet != "Baseline"))


# Diet model ####
lm(breast_cancer ~ diet + 0, data = mouse) -> model
model |> summary() |> coef() |> as.data.frame() -> c
rownames(c) |> str_remove("diet") -> c$diet

qt(0.975, df = df.residual(model)) -> k
c |> filter(diet != "Baseline") |> ggplot(
  aes(
    x = diet,
    fill = diet,
    y = Estimate,
    ymin = Estimate - k * `Std. Error`,
    ymax = Estimate + k * `Std. Error`
  )
) +
  ylim(-0.05, 0.55) -> p

print(
  p + geom_col() + geom_errorbar() +
    ggtitle("Rate of Breast Cancer for each Diet",
            subtitle = "95% confidence intervals")
)

# Combination of age and diet model ####
age <- c(
  "49" = "49/50",
  "50" = "49/50",
  "81" = "81/82",
  "82" = "81/82",
  "17" = "17/18",
  "18" = "17/18",
  "10" = "10"
)
mouse$age <- age[mouse$sacrifice_week |> as.character()]

lm(breast_cancer ~ age:diet + 0,
   data = mouse,
   subset = sacrifice_week > 40) -> model0
model0 |> summary() |> coef() |> as.data.frame() -> c0
rownames(c0) |> str_remove("diet") |> str_remove("age") -> c0$diet

qt(0.975, df = df.residual(model0)) -> k0
c0 |> filter(diet != "Baseline") |> ggplot(
  aes(
    x = diet,
    fill = diet,
    y = Estimate,
    ymin = Estimate - k0 * `Std. Error`,
    ymax = Estimate + k0 * `Std. Error`
  )
) +
  ylim(-0.05, 0.55) -> p0

print(
  p0 + geom_col() + geom_errorbar() +
    ggtitle("Rate of Breast Cancer for each Diet",
            subtitle = "95% confidence intervals")
)


# merge 'ICRR' and 'ICRRF' ####

mouse$d <- ifelse(mouse$diet == "ICR-R" | mouse$diet == "ICR-RF",
                  "ICR", mouse$diet)
features$diet <- as.factor(features$diet)
features$diet <- relevel(features$diet, ref = "BASELINE")

lm(breast_cancer ~ d + 0, data = mouse) -> model1
model1 |> summary() |> coef() |> as.data.frame() -> c1
rownames(c1) |> str_remove("d") -> c1$d
qt(0.975, df = df.residual(model1)) -> k1
c1 |>
  filter(diet != "Baseline") |>
  ggplot(
    aes(
      x = diet,
      fill = diet,
      y = Estimate,
      ymin = Estimate - k1 * `Std. Error`,
      ymax = Estimate + k1 * `Std. Error`
    )
  ) +
  ylim(-0.05, 0.55) -> p1

print(
  p1 + geom_col() + geom_errorbar(width = 0.2) +
    ggtitle("Rate of Breast Cancer for each Diet",
            subtitle = "95% confidence intervals")
)


mouse |> filter(sacrifice_week > 0) -> mouse0

mouse0$diet <- as.factor(mouse0$diet)
mouse0$diet <- relevel(mouse0$diet, ref = "Baseline")

lm(breast_cancer ~ d, data = mouse) -> model01
