## code to prepare `respiratory` dataset goes here

library(dplyr)
library(rio)

path <- file.path("data-raw", "respiratory.csv")

respiratory <- rio::import(path) |>
  transmute(
    id = as.integer(rep(1:111, each = 4)),
    visit = as.integer(visit),
    status = as.integer(outcome),
    treatment = factor(
      if_else(treat == "P", "placebo", "active"),
      levels = c("placebo", "active")
    ),
    baseline = as.integer(baseline),
    gender = factor(
      if_else(sex == "M", "male", "female"),
      levels = c("female", "male")
    ),
    age = as.numeric(age),
    center = factor(
      if_else(center == 1, "C1", "C2"),
      levels = c("C1", "C2")
    )
  )

stopifnot(
  nrow(respiratory) == 444L,
  all(respiratory$status %in% c(0L, 1L)),
  all(respiratory$baseline %in% c(0L, 1L)),
  !anyNA(respiratory$treatment),
  !anyNA(respiratory$gender),
  !anyNA(respiratory$center)
)

usethis::use_data(respiratory, overwrite = TRUE, compress = "xz")
