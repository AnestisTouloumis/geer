## code to prepare `epilepsy` dataset goes here

library("ALA")
library("tidyverse")
epilepsy <-
  ALA::epilepsy |>
  mutate(response = if_else(nSeizures <= 5, 1, 0)) |>
  relocate(id, week, response, treatment, age, nSeizures) |>
  rename(y = response) |>
  mutate(id = factor(id),
         week = as.numeric(week),
         y = as.numeric(y),
         treatment = as.factor(treatment),
         age = as.numeric(age),
         nSeizures = as.numeric(nSeizures)) |>
  as.data.frame()

rownames(epilepsy) <- 1:nrow(epilepsy)
usethis::use_data(epilepsy, overwrite = TRUE)
