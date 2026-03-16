## code to prepare `depression` dataset goes here

library(dplyr)
library(rio)

path <- file.path("data-raw", "depression.csv")

raw <- rio::import(path)

stopifnot(all(raw$group %in% c("placebo", "oestrogen")))

baseline_tbl <- raw |>
  filter(visit == -1) |>
  transmute(
    subj = as.integer(subj),
    baseline = as.numeric(dep)
  )

stopifnot(
  all(table(baseline_tbl$subj) == 1L),
  !anyNA(baseline_tbl$baseline)
)

depression <- raw |>
  filter(visit != -1) |>
  mutate(
    subj = as.integer(subj),
    visit = as.integer(visit)
  ) |>
  left_join(baseline_tbl, by = "subj") |>
  transmute(
    id = subj,
    visit = visit,
    score = as.numeric(dep),
    treatment = factor(group, levels = c("placebo", "oestrogen")),
    baseline = as.numeric(baseline)
  )

stopifnot(
  nrow(depression) == 366L,
  !anyNA(depression$baseline),
  !anyNA(depression$treatment)
)

usethis::use_data(depression, overwrite = TRUE, compress = "xz")
