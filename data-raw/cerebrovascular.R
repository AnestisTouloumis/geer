## code to prepare `cerebrovascular` dataset goes here

library("tidyverse")
library("rio")
cerebrovascular <-
  rio::import("./data-raw/cerebrovascular.csv") |>
  mutate(treatment =
           if_else(sequence == "P->A" & period == 1,
                   "placebo",
                   if_else(sequence == "P->A" & period == 2,
                           "active",
                           if_else(sequence == "A->P" & period == 1,
                                   "active",
                                   "placebo")
                   )
           ),
         ecg = if_else(ecg == "normal", 1, 0)
  ) |>
  select(-sequence) |>
  mutate(id = as.numeric(id),
         period = as.numeric(period),
         ecg = as.numeric(ecg),
         treatment = factor(treatment)) |>
  relocate(id, period, ecg, treatment)
rownames(cerebrovascular) <- 1:nrow(cerebrovascular)
cerebrovascular <- as_tibble(cerebrovascular)
usethis::use_data(cerebrovascular, overwrite = TRUE)
