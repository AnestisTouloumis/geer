## code to prepare `epilepsy` dataset goes here

library("ALA")
library("tidyverse")
epilepsy <-
  epilepsy |>
  mutate(response = if_else(nSeizures <= 5, 1, 0))
