#!/usr/bin/env Rscript
library(limma)
library(stringr)

## TODO: drop Brain

# input
MOUSE_DATA <- paste(RESULT_FOLDER, "mouse-data-features.RData", sep = "/")

# output
FITTED_FILE <- "fit.eb.RData"

# read `mouse.Data`, `features`, `mouse_4.1`
load(MOUSE_DATA)

# Evaluate absolute expression for each diet:age:tissue
# design <- model.matrix(~ 0 + diet:age:tissue, data = features)

# TODO: compare diet depending on tissue independent of age
# Evaluate absolute expression for each diet:tissue
design <- model.matrix(~ 0 + diet:age:tissue, data = features)
design <- design[, colSums(design)>0]
colnames(design) |> 
  str_remove("diet") |>
  str_remove("tissue") |>
  str_remove("age") |> 
  str_replace_all("/", ".") |> 
  str_replace_all(":", ".") -> colnames(design)
# now we have `design`

# we fit our normalized expression data to our design
fit <- lmFit(mouse.Data, design)

## now we can do the contrasts that we choose
## we compare each condition to the baseline
## at baseline (day 0) there is no diet application
## in the first age group, for 50 weeks diet is applied (from day 0)
## in the second age group, for 81 weeks diet is applied (from day 0)
## finally we decided to merge two age groups
## basically we ignore the effect of age in gene expression
contrasts <- makeContrasts(levels = design,
  ## Blood
  d.AL.49.50.Blood  = AL.49.50.Blood  - BASELINE.10.Blood,
  d.CCR.49.50.Blood = CCR.49.50.Blood - BASELINE.10.Blood,
  # d.CCR.49.50.Blood = CCR.49.50.Blood - AL.49.50.Blood  - BASELINE.10.Blood
  d.ICR.49.50.Blood = ICR.49.50.Blood - BASELINE.10.Blood,
  d.AL.81.82.Blood  = AL.81.82.Blood  - BASELINE.10.Blood,
  d.CCR.81.82.Blood = CCR.81.82.Blood - BASELINE.10.Blood,
  d.ICR.81.82.Blood = ICR.81.82.Blood - BASELINE.10.Blood,
  ## Brain
  d.AL.49.50.Brain  = AL.49.50.Brain  - BASELINE.10.Brain,
  d.CCR.49.50.Brain = CCR.49.50.Brain - BASELINE.10.Brain,
  d.ICR.49.50.Brain = ICR.49.50.Brain - BASELINE.10.Brain,
  d.AL.81.82.Brain  = AL.81.82.Brain  - BASELINE.10.Brain,
  d.CCR.81.82.Brain = CCR.81.82.Brain - BASELINE.10.Brain,
  d.ICR.81.82.Brain = ICR.81.82.Brain - BASELINE.10.Brain,
  ## MFP
  d.AL.49.50.MFP    = AL.49.50.MFP    - BASELINE.10.MFP,
  d.CCR.49.50.MFP   = CCR.49.50.MFP   - BASELINE.10.MFP,
  d.ICR.49.50.MFP   = ICR.49.50.MFP   - BASELINE.10.MFP,
  d.AL.81.82.MFP    = AL.81.82.MFP    - BASELINE.10.MFP,
  d.CCR.81.82.MFP   = CCR.81.82.MFP   - BASELINE.10.MFP,
  d.ICR.81.82.MFP   = ICR.81.82.MFP   - BASELINE.10.MFP
  # d.CR = CCR.49.50.Blood - AL.49.50.Blood + BASELINE.10.Blood
  # d.IR = ICR.49.50.Blood - AL.49.50.Blood + BASELINE.10.Blood
  # d.RR = ICR.49.50.Blood - CCR.49.50.Blood + BASELINE.10.Blood
)

fit.cont <- contrasts.fit(fit, contrasts)
fit.eb <- eBayes(fit.cont)

# `fit.eb` contains all the results we care
save(fit.eb, file = FITTED_FILE)
