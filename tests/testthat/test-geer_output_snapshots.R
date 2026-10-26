list2env(get_cached_test_fixtures("geer_output_snapshots"), environment())

skip_on_cran()




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
