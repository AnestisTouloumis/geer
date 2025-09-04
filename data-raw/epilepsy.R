## code to prepare `epilepsy` dataset goes here

library("rio")
library("tidyverse")
epilepsy <-
  rio::import("./data-raw/epilepsy.csv") |>
  mutate(id = as.numeric(id),
         visit = as.numeric(time),
         seizures = as.numeric(seizures),
         treatment = if_else(treatment == "placebo", "placebo", "progabide") |>
                     as.factor(),
         lnage = log(as.numeric(age)),
         lnbaseline = log(as.numeric(base)/4)) |>
  select(id, visit, seizures, treatment, lnbaseline, lnage)
rownames(epilepsy) <- 1:nrow(epilepsy)
epilepsy <- as_tibble(epilepsy)
usethis::use_data(epilepsy, overwrite = TRUE)
