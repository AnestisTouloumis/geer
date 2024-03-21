## code to prepare `leprosy` dataset goes here

library("tidyverse")
library("ALA")
data("leprosy")
leprosy <-
  ALA::leprosy |>
  rename(bacilli = nBacilli,
         treatment = drug) |>
  relocate(id, period, bacilli, treatment)

rownames(leprosy) <- 1:nrow(leprosy)

usethis::use_data(leprosy, overwrite = TRUE)
