## code to prepare `rinse` dataset goes here

library(dplyr)
library(rio)

path <- file.path("data-raw", "rinse.csv")

raw <- rio::import(path)

baseline_tbl <- raw |>
  filter(time == 0) |>
  transmute(
    subject = as.integer(subject),
    baseline = as.numeric(score)
  )

stopifnot(
  all(table(baseline_tbl$subject) == 1L),
  !anyNA(baseline_tbl$baseline)
)

rinse <- raw |>
  filter(time != 0) |>
  mutate(
    subject = as.integer(subject),
    time = as.integer(time)
  ) |>
  left_join(baseline_tbl, by = "subject") |>
  transmute(
    id = subject,
    time = time,
    score = as.numeric(score),
    treatment = factor(
      case_when(
        rinse == "Placebo" ~ "placebo",
        TRUE ~ as.character(rinse)
      ),
      levels = c("placebo", "A", "B")
    ),
    baseline = as.numeric(baseline),
    gender = factor(
      if_else(gender == "Male", "male", "female"),
      levels = c("female", "male")
    ),
    age = as.numeric(age),
    smoke = factor(
      if_else(smoke == "Yes", "yes", "no"),
      levels = c("no", "yes")
    )
  )

stopifnot(
  nrow(rinse) == 218L,
  !anyNA(rinse$treatment),
  !anyNA(rinse$baseline),
  !anyNA(rinse$gender),
  !anyNA(rinse$smoke)
)

usethis::use_data(rinse, overwrite = TRUE, compress = "xz")
