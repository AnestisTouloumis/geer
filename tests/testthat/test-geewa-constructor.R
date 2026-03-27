testthat::local_edition(3)

fit <- geewa(
  seizures ~ treatment + lnbaseline + lnage,
  data = test_data$epilepsy,
  id = id,
  repeated = visit
)

test_that("geer.default reconstructs a valid geer object", {
  out <- geer.default(unclass(fit))

  expect_s3_class(out, "geer")
  expect_equal(coef(out), coef(fit))
})

test_that("new_geer validates required components and fills optional metadata", {
  x_missing <- unclass(fit)
  x_missing$coefficients <- NULL

  expect_error(
    new_geer(x_missing),
    "missing required components"
  )

  x_optional <- unclass(fit)
  x_optional <- x_optional[setdiff(names(x_optional), c("na.action", "contrasts", "xlevels"))]

  out <- new_geer(x_optional)

  expect_s3_class(out, "geer")
  expect_null(out$na.action)
  expect_null(out$contrasts)
  expect_equal(out$xlevels, list())
})

test_that("validate_geer accepts a valid geer object", {
  expect_s3_class(validate_geer(fit), "geer")
})

test_that("validate_geer checks key structural consistency", {
  bad_cov <- new_geer(unclass(fit))
  bad_cov$naive_covariance <- matrix(0, 2, 2)

  expect_error(
    validate_geer(bad_cov),
    "naive_covariance"
  )

  bad_clusters <- new_geer(unclass(fit))
  bad_clusters$clusters_no <- bad_clusters$clusters_no + 1

  expect_error(
    validate_geer(bad_clusters),
    "clusters_no"
  )

  bad_coef <- fit
  names(bad_coef$coefficients) <- NULL

  expect_error(
    validate_geer(bad_coef),
    "'coefficients' must be a named numeric vector"
  )

  bad_x <- fit
  bad_x$x <- bad_x$x[, -1, drop = FALSE]

  expect_error(
    validate_geer(bad_x),
    "number of columns of 'x' must match length of 'coefficients'"
  )

  bad_dimnames <- fit
  bad_dimnames$robust_covariance <- diag(length(bad_dimnames$coefficients))
  rownames(bad_dimnames$robust_covariance) <- names(bad_dimnames$coefficients)
  colnames(bad_dimnames$robust_covariance) <- paste0("b", seq_along(bad_dimnames$coefficients))

  expect_error(
    validate_geer(bad_dimnames),
    "dimnames of 'robust_covariance' must match names of 'coefficients'"
  )
})
