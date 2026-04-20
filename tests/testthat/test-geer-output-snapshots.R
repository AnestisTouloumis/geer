testthat::local_edition(3)

data(epilepsy, package = "geer")

fit_full <- geewa(
  seizures ~ treatment + lnbaseline + lnage,
  data = epilepsy,
  id = id,
  repeated = visit
)

fit_reduced <- geewa(
  seizures ~ treatment,
  data = epilepsy,
  id = id,
  repeated = visit
)

test_that("core printed outputs remain stable", {
  expect_snapshot_output(print(fit_full))
  expect_snapshot_output(print(summary(fit_full)))
  expect_snapshot_output(
    print(anova(fit_full, test = "wald", cov_type = "robust"))
  )
})


test_that("stepwise printed outputs remain stable", {
  expect_snapshot_output(
    print(add1(
      fit_reduced,
      scope = ~ treatment + lnbaseline + lnage,
      test = "wald",
      cov_type = "robust"
    ))
  )
  expect_snapshot_output(
    print(drop1(fit_full, test = "wald", cov_type = "robust"))
  )
})
