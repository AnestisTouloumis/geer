## code to prepare `cerebrovascular` dataset goes here

library("tidyverse")
library("ALA")
cerebrovascular <-
  ALA::ecg |>
  mutate(treatment =
           if_else(sequence == "P->A" & period == 1,
                   "placebo",
                   if_else(sequence == "P->A" & period == 2,
                           "active",
                           if_else(sequence == "A->P" & period == 1,
                                   "active",
                                   "placebo")
                   )
           )
  ) |>
  select(-sequence) |>
  relocate(id, period, ecg, treatment) |>
  mutate(id = factor(id),
         period = factor(period),
         ecg = factor(ecg),
         treatment = factor(treatment)) |>
  as.data.frame()
rownames(cerebrovascular) <- 1:nrow(cerebrovascular)
usethis::use_data(cerebrovascular, overwrite = TRUE)
