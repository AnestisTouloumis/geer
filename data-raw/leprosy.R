## code to prepare `leprosy` dataset goes here

library("tidyverse")
library("ALA")
data("leprosy")
leprosy <-
  ALA::leprosy |>
  rename(bacilli = nBacilli,
         treatment = drug) |>
  mutate(time = if_else(period == "post", 1, 0)) |>
  relocate(id, period, time, bacilli, treatment)
rownames(leprosy) <- 1:nrow(leprosy)

usethis::use_data(leprosy, overwrite = TRUE)


fitted_model_gee <- geewa_binary(
  formula = I(obesity == "Obese") ~ age + zygosity + ancestry + bacteroides,
