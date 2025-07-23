## code to prepare `epilepsy` dataset goes here

library("rio")
library("tidyverse")
epilepsy <-
  rio::import("./data-raw/seizures.csv") |>
  mutate(id = as.numeric(id),
         visit = as.numeric(time),
         seizures = as.numeric(seizures),
         treatment = if_else(treatment == "Placebo", "placebo", "progabide") |>
                     as.factor(),
         age = as.numeric(age),
         base = as.numeric(base)) |>
  select(id, visit, seizures, treatment, base, age)
rownames(epilepsy) <- 1:nrow(epilepsy)
epilepsy <- as_tibble(epilepsy)
usethis::use_data(epilepsy, overwrite = TRUE)
