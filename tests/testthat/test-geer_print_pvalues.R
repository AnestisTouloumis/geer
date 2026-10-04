test_that("fix_small_pvalue_labels relabels tiny p-values and keeps width", {
  lines <- c("a  1.0  <1e-04 ***", "b  2.0  < 1e-04 ***", "c  3.0  0.0123 *")
  out <- fix_small_pvalue_labels(lines)
  expect_equal(out[3L], lines[3L])
  expect_true(all(grepl("<0.0001", out[1:2], fixed = TRUE)))
  expect_false(any(grepl("1e-04", out, fixed = TRUE)))
  expect_equal(nchar(out[2L]), nchar(lines[2L]))
})

test_that("summary and anova print p-values below 0.0001 as <0.0001", {
  coefs <- cbind(
    Estimate = c(1, 2), `Std. Error` = c(0.1, 0.1),
    `z value` = c(10, 20), `Pr(>|z|)` = c(1e-20, 0.2)
  )
  rownames(coefs) <- c("x1", "x2")
  out <- utils::capture.output(print_with_small_pvalues(function() {
    stats::printCoefmat(coefs, eps.Pvalue = 1e-04)
  }))
  expect_true(any(grepl("<0.0001", out, fixed = TRUE)))
  expect_false(any(grepl("e-", out, fixed = TRUE)))

  tab <- structure(
    data.frame(Df = 1, Chi = 90, `Pr(>Chi)` = 1e-15, check.names = FALSE),
    heading = "Test", class = c("geer_anova", "anova", "data.frame")
  )
  out <- utils::capture.output(res <- print(tab))
  expect_true(any(grepl("<0.0001", out, fixed = TRUE)))
  expect_identical(res, tab)
})

test_that("htest results print p-values below 0.0001 as < 0.0001", {
  small <- structure(
    list(statistic = c(X = 90), parameter = c(df = 1), p.value = 3e-21,
         method = "Demo test", data.name = "x"),
    class = c("geer_htest", "htest")
  )
  out <- utils::capture.output(res <- print(small))
  expect_true(any(grepl("p-value < 0.0001", out, fixed = TRUE)))
  expect_false(any(grepl("e-", out, fixed = TRUE)))
  expect_identical(res, small)
  moderate <- small
  moderate$p.value <- 0.0123
  out <- utils::capture.output(print(moderate))
  expect_true(any(grepl("p-value = 0.0123", out, fixed = TRUE)))
  boundary <- small
  boundary$p.value <- 1e-04
  out <- utils::capture.output(print(boundary))
  expect_false(any(grepl("< 0.0001", out, fixed = TRUE)))
})

test_that("the runs and MCAR tests carry the geer_htest class", {
  expect_s3_class(runs_test(fit_geewa_pois_exch), "geer_htest")
})
