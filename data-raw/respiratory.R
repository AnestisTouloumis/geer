## code to prepare `respiratory` dataset goes here

library("rio")
library("tidyverse")
respiratory <-
  rio::import("./data-raw/respiratory.csv") |>
  mutate(treat = if_else(treat == "P", "placebo", "active"),
         sex = if_else(sex == "M", "male", "female"),
         center = if_else(center == 1, "C1", "C2"),
         id = rep(1:111, each = 4)) |>
  mutate(id = as.numeric(id),
         visit = as.numeric(visit),
         y = as.numeric(outcome),
         treatment = as.factor(treat),
         baseline = as.numeric(baseline),
         age = as.numeric(age),
         gender = factor(sex),
         center = factor(center)) |>
  select(id, visit, y, treatment, baseline, gender, age, center)
rownames(respiratory) <- 1:nrow(respiratory)
respiratory <- as_tibble(respiratory)
usethis::use_data(respiratory, overwrite = TRUE)
