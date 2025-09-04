## code to prepare `leprosy` dataset goes here

library("tidyverse")
library("rio")
leprosy <-
  rio::import("./data-raw/leprosy.csv") |>
  rename(bacilli = nBacilli,
         treatment = drug) |>
  mutate(id = as.numeric(id),
         period = factor(period),
         bacilli = as.numeric(bacilli),
         treatment = factor(treatment)) |>
  relocate(id, period, bacilli, treatment)
rownames(leprosy) <- 1:nrow(leprosy)
leprosy <- as_tibble(leprosy)
usethis::use_data(leprosy, overwrite = TRUE)
