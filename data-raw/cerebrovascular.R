## code to prepare `cerebrovascular` dataset goes here

library(dplyr)
library(rio)

path <- file.path("data-raw", "cerebrovascular.csv")

cerebrovascular <- rio::import(path) |>
  transmute(
    id = as.integer(id),
    period = as.integer(period),
    ecg = as.integer(if_else(ecg == "normal", 1L, 0L)),
    treatment = factor(
      case_when(
        sequence == "P->A" & period == 1 ~ "placebo",
        sequence == "P->A" & period == 2 ~ "active",
        sequence == "A->P" & period == 1 ~ "active",
        sequence == "A->P" & period == 2 ~ "placebo",
        TRUE ~ NA_character_
      ),
      levels = c("placebo", "active")
    )
  )

stopifnot(
  all(c("id", "period", "ecg", "treatment") %in% names(cerebrovascular)),
  all(cerebrovascular$ecg %in% c(0L, 1L)),
  !anyNA(cerebrovascular$treatment)
)

usethis::use_data(cerebrovascular, overwrite = TRUE, compress = "xz")
