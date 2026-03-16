## code to prepare `leprosy` dataset goes here

library(dplyr)
library(rio)

path <- file.path("data-raw", "leprosy.csv")

leprosy <- rio::import(path) |>
  transmute(
    id = as.integer(id),
    period = factor(period, levels = c("pre", "post")),
    bacilli = as.integer(nBacilli),
    treatment = factor(drug, levels = c("A", "B", "C"))
  )

stopifnot(
  all(c("id", "period", "bacilli", "treatment") %in% names(leprosy)),
  !anyNA(leprosy$period),
  !anyNA(leprosy$treatment)
)

usethis::use_data(leprosy, overwrite = TRUE, compress = "xz")
