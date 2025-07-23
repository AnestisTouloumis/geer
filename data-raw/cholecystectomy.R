## code to prepare `cholecystectomy` dataset goes here

library("tidyverse")
library("rio")
cholecystectomy <-
  rio::import("./data-raw/cholecystectomy.csv") |>
  select(-pain) |>
  mutate(treatment = if_else(treatment == "A", "active", "placebo"),
         gender = if_else(gender == "F", "female", "male"),
         treatment = factor(treatment),
         gender = factor(gender),
         age = as.numeric(age),
         time = as.numeric(time),
         pain2 = as.numeric(pain2)) |>
  rename(pain = pain2) |>
  select(id, time, pain, treatment, gender, age)
rownames(cholecystectomy) <- 1:nrow(cholecystectomy)
cholecystectomy <- as_tibble(cholecystectomy)


usethis::use_data(cholecystectomy, overwrite = TRUE)
