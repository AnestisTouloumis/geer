## code to prepare `depression` dataset goes here

library("rio")
library("tidyverse")
depression <-
  rio::import("./data-raw/depression.csv") |>
  group_by(subj) |>
  mutate(base = dep[visit==-1]) |>
  ungroup() |>
  filter(visit != -1) |>
  mutate(id = as.numeric(subj),
         visit = as.numeric(visit),
         treatment = factor(group),
         score = as.numeric(dep),
         baseline = as.numeric(base)) |>
  select(id, visit, score, treatment, baseline)
depression <- as_tibble(depression)
usethis::use_data(depression, overwrite = TRUE)
