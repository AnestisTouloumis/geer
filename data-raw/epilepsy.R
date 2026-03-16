## code to prepare `epilepsy` dataset goes here

library(dplyr)
library(rio)

path <- file.path("data-raw", "epilepsy.csv")

epilepsy <- rio::import(path) |>
  transmute(
    id = as.integer(id),
    visit = as.integer(time),
    seizures = as.integer(seizures),
    treatment = factor(
      case_when(
        treatment == "placebo" ~ "placebo",
        treatment == "Progabide" ~ "progabide"
      )
      ),
    lnbaseline = log(as.numeric(base) / 4),
    lnage = log(as.numeric(age))
  )

stopifnot(
  !anyNA(epilepsy),
  all(is.finite(epilepsy$lnbaseline)),
  all(is.finite(epilepsy$lnage))
)

usethis::use_data(epilepsy, overwrite = TRUE, compress = "xz")
