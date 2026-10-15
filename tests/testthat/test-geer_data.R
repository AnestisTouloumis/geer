test_that("package datasets load and contain essential columns", {
  expected_columns <- list(
    cerebrovascular = c("id", "period", "ecg", "treatment"),
    cholecystectomy = c("id", "time", "pain", "treatment", "gender", "age"),
    depression = c("id", "visit", "score", "treatment", "baseline"),
    epilepsy = c("id", "visit", "seizures", "treatment", "lnbaseline", "lnage"),
    leprosy = c("id", "period", "bacilli", "treatment"),
    respiratory = c("id", "visit", "status", "treatment"),
    rinse = c("id", "time", "score", "treatment", "baseline")
  )
  for (nm in names(expected_columns)) {
    dat <- load_geer_data(nm)
    expect_s3_class(dat, "data.frame")
    expect_true(all(expected_columns[[nm]] %in% names(dat)))
    expect_gt(nrow(dat), 0L)
  }
})
