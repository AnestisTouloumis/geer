## code to prepare `epilepsy` dataset goes here

library("glmtoolbox")
library("tidyverse")
data("Seizures")
epilepsy <-
  Seizures |>
  relocate(id, time, seizures, treatment, base, age) |>
  mutate(id = factor(id),
         time = as.numeric(time),
         seizures = as.numeric(seizures),
         treatment = as.factor(treatment),
         age = log(as.numeric(age)),
         base = log(as.numeric(base)/4),
         visit4 = if_else(time == 4, 1, 0)) |>
  as.data.frame()

rownames(epilepsy) <- 1:nrow(epilepsy)
usethis::use_data(epilepsy, overwrite = TRUE)


select(subject, age, treatment, base, period, seizure.rate) |>
  mutate(base = log(base),
         age = log(age),
         visit4 = I(period == 4),
         per = -log(4)) |>
  select(subject, age, treatment, base, visit4, period, seizure.rate, per)
