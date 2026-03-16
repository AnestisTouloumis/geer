## code to prepare `cholecystectomy` dataset goes here

library(dplyr)
library(rio)

path <- file.path("data-raw", "cholecystectomy.csv")

cholecystectomy <- rio::import(path) |>
  transmute(
    id = as.integer(id),
    time = as.integer(time),
    pain = as.integer(pain2),  # 0/1 indicator, as documented
    treatment = factor(
      if_else(treatment == "A", "active", "placebo"),
      levels = c("placebo", "active")
    ),
    gender = factor(
      if_else(gender == "F", "female", "male"),
      levels = c("female", "male")
    ),
    age = as.numeric(age)
  )

stopifnot(
  nrow(cholecystectomy) == 246L,
  all(cholecystectomy$pain %in% c(0L, 1L)),
  !anyNA(cholecystectomy$treatment),
  !anyNA(cholecystectomy$gender)
)

usethis::use_data(cholecystectomy, overwrite = TRUE, compress = "xz")
