## code to prepare `cerebrovascular` dataset goes here

usethis::use_data(cerebrovascular, overwrite = TRUE)
library("tidyverse")
library("ALA")
data("ecg")
cerebrovascular <-
  ecg |>
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
  select(-sequence)
