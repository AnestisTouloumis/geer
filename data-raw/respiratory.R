## code to prepare `respiratory` dataset goes here

library("tidyverse")
data("respiratory", package = "mdhglm")
respiratory <-
  respiratory |>
  mutate(id = patient) |>
  group_by(id) |>
  mutate(visit = 1:4,
         gender = sex) |>
  ungroup(id) |>
  select(-patient, -past, - sex) |>
  relocate(id, y, baseline, treatment, visit, gender, age, center)
usethis::use_data(respiratory, overwrite = TRUE)
