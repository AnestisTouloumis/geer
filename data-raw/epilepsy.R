## code to prepare `epilepsy` dataset goes here

library("rio")
library("tidyverse")
epilepsy <-
  rio::import("./data-raw/Seizures.csv") |>
  relocate(id, time, seizures, treatment, base, age) |>
  mutate(id = factor(id),
         visit = as.numeric(time),
         seizures = as.numeric(seizures),
         treatment = as.factor(treatment),
         age = age,
         base = base) |>
  select(id, visit, seizures, treatment, base, age) |>
  as.data.frame()

rownames(epilepsy) <- 1:nrow(epilepsy)
usethis::use_data(epilepsy, overwrite = TRUE)
