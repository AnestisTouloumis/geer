## code to prepare `rinse` dataset goes here

library("rio")
library("tidyverse")
rinse <-
  rio::import("./data-raw/rinse.csv") |>
  group_by(subject) |>
  mutate(base = score[time == 0]) |>
  ungroup() |>
  filter(time != 0) |>
  mutate(rinse = if_else(rinse == "Placebo", "placebo", rinse),
         gender = if_else(gender == "Male", "male", "female"),
         smoke = if_else(smoke == "Yes", "yes", "no")
         ) |>
  mutate(id = as.numeric(subject),
         time = as.numeric(time),
         score = as.numeric(score),
         treatment = as.factor(rinse),
         gender = as.factor(gender),
         baseline = as.numeric(base),
         age = as.numeric(age),
         smoke = as.factor(smoke)) |>
  select(id, time, score, treatment, baseline, gender, age, smoke)
rinse <- as_tibble(rinse)
usethis::use_data(rinse, overwrite = TRUE)
