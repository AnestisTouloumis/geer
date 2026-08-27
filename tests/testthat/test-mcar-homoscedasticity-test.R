testthat::local_edition(3)


make_jj_screening_data <- function() {
  set.seed(9104)
  x <- matrix(stats::rnorm(180), nrow = 60, ncol = 3)
  colnames(x) <- c("visit1", "visit2", "visit3")
  x[31:45, 3] <- NA_real_
  x[46:60, 2:3] <- NA_real_
  x
}


test_that("internal Anderson-Darling calculation follows Scholz-Stephens", {
  x <- c(
    -1.2, -0.5, 0.1, 0.4, 0.9, 1.1, 1.5,
    -1.0, -0.2, 0.05, 0.6, 0.8, 1.3, 1.8, 2.1,
    -0.8, -0.1, 0.2, 0.45, 0.7, 1.0, 1.4, 1.7, 2.0
  )
  group <- rep(1:3, c(7, 8, 9))
  out <- geer:::jj_anderson_darling_test(x, group, c(7L, 8L, 9L))

  expect_equal(out$statistic, 0.8000654142, tolerance = 1e-9)
  expect_equal(
    out$group.statistics,
    c(0.40675810, 0.20711942, 0.18618789),
    tolerance = 1e-7
  )
  expect_equal(out$variance, 0.9526947675, tolerance = 1e-9)
  expect_equal(out$standardized, -1.229364538, tolerance = 1e-8)
  expect_true(is.finite(out$p.value))
  expect_true(out$p.value > 0 && out$p.value < 1)
})


test_that("internal modified Hawkins calculation follows the published transformation", {
  g1 <- cbind(
    seq(-1, 1, length.out = 8),
    seq(-0.8, 1.2, length.out = 8) +
      c(0, 0.1, -0.1, 0.05, -0.05, 0.08, -0.08, 0)
  )
  g2 <- cbind(
    seq(-0.7, 1.3, length.out = 8),
    seq(-1.1, 0.9, length.out = 8) +
      c(0.05, -0.1, 0.1, -0.02, 0.07, -0.04, 0.09, -0.06)
  )
  g3 <- cbind(
    seq(-1.2, 0.8, length.out = 8),
    seq(-0.9, 1.1, length.out = 8) +
      c(-0.07, 0.04, 0.08, -0.09, 0.03, 0.1, -0.05, 0.02)
  )
  completed <- rbind(g1, g2, g3)
  group <- rep(1:3, each = 8)

  out <- geer:::jj_hawkins_test(
    completed,
    group = group,
    group_counts = c(8L, 8L, 8L),
    nrep = 20L,
    n_min = 2L
  )

  expect_equal(
    out$pooled.covariance,
    matrix(c(
      0.48979592, 0.48809524,
      0.48809524, 0.49173920
    ), 2, 2),
    tolerance = 1e-7
  )
  expect_equal(
    out$group.statistics,
    c(5.705201788, 3.259474395, 5.668594317),
    tolerance = 1e-8
  )
  expect_equal(out$statistic, 7.314033418, tolerance = 1e-8)
  expect_identical(out$parameter, 6L)
  expect_equal(out$p.value, 0.2927792548, tolerance = 1e-9)
})


test_that("mcar_homoscedasticity_test returns both diagnostics in auto mode", {
  x <- make_jj_screening_data()
  set.seed(212)
  rng_before <- .Random.seed

  out <- mcar_homoscedasticity_test(
    x,
    method = "auto",
    n_min = 2L,
    seed = 110L
  )

  expect_s3_class(out, "htest")
  expect_s3_class(out, "mcar_homoscedasticity_test")
  expect_identical(out$method.requested, "auto")
  expect_identical(out$imputation.requested, "distribution-free")
  expect_identical(out$imputation.used, "distribution-free")
  expect_equal(out$complete.cases, 30L)
  expect_equal(out$n, 60L)
  expect_equal(out$p, 3L)
  expect_equal(nrow(out$tests), 2L)
  expect_identical(out$tests$test, c("hawkins", "nonparametric"))
  expect_true(all(is.finite(out$tests$statistic)))
  expect_true(all(out$tests$p.value >= 0 & out$tests$p.value <= 1))
  expect_equal(out$pattern.counts, c(pattern1 = 30L, pattern2 = 15L, pattern3 = 15L))
  expect_equal(dim(out$patterns), c(3L, 3L))
  expect_equal(
    unname(out$patterns),
    matrix(c(
      0L, 0L, 0L,
      0L, 0L, 1L,
      0L, 1L, 1L
    ), nrow = 3L, byrow = TRUE)
  )
  expect_false(anyNA(out$imputed.data))
  expect_true(out$selected.test %in% c("hawkins", "nonparametric"))
  expect_identical(.Random.seed, rng_before)
})


test_that("nonparametric mode does not run the Neyman uniformity test", {
  x <- make_jj_screening_data()
  out <- mcar_homoscedasticity_test(
    x,
    method = "nonparametric",
    seed = 123L
  )

  expect_identical(out$selected.test, "nonparametric")
  expect_identical(out$tests$test, "nonparametric")
  expect_true(is.na(out$hawkins$p.value))
  expect_true(all(is.na(out$hawkins$group.p.values)))
  expect_equal(out$p.value, out$nonparametric$p.value)
})


test_that("distribution-free imputation falls back to normal imputation when needed", {
  set.seed(490)
  x <- matrix(stats::rnorm(72), nrow = 24, ncol = 3)
  x[9:16, 3] <- NA_real_
  x[17:24, 2:3] <- NA_real_

  expect_warning(
    out <- mcar_homoscedasticity_test(
      x,
      method = "nonparametric",
      seed = 10L
    ),
    "normal-theory imputation is used"
  )
  expect_identical(out$imputation.requested, "distribution-free")
  expect_identical(out$imputation.used, "normal")
  expect_equal(out$complete.cases, 8L)
})


test_that("small missingness patterns are omitted", {
  set.seed(901)
  x <- matrix(stats::rnorm(120), nrow = 40, ncol = 3)
  x[21:30, 3] <- NA_real_
  x[31:37, 2:3] <- NA_real_
  x[38:40, 1] <- NA_real_

  expect_warning(
    out <- mcar_homoscedasticity_test(
      x,
      method = "nonparametric",
      imputation = "normal",
      seed = 9L
    ),
    "fewer than 7 cases"
  )

  expect_equal(out$n, 37L)
  expect_equal(nrow(out$omitted.patterns), 1L)
  expect_equal(out$omitted.patterns$n, 3L)
})


test_that("mcar_homoscedasticity_test reconstructs responses from a geer fit", {
  set.seed(5502)
  subjects <- 45L
  visits <- 3L
  id <- rep(seq_len(subjects), each = visits)
  visit <- rep(seq_len(visits), times = subjects)
  trt <- rep(rep(c(0, 1), length.out = subjects), each = visits)
  y <- 1 + 0.3 * trt + 0.2 * visit + stats::rnorm(length(id))

  y[id %in% 26:35 & visit == 3] <- NA_real_
  y[id %in% 36:45 & visit %in% 2:3] <- NA_real_
  dat <- data.frame(id = id, visit = visit, trt = trt, y = y)

  fit <- geewa(
    y ~ trt + factor(visit),
    data = dat,
    id = id,
    repeated = visit,
    family = gaussian(),
    corstr = "independence"
  )

  wide <- matrix(NA_real_, nrow = subjects, ncol = visits)
  wide[cbind(dat$id, dat$visit)] <- dat$y

  from_fit <- mcar_homoscedasticity_test(
    fit,
    method = "nonparametric",
    imputation = "normal",
    seed = 77L
  )
  from_wide <- mcar_homoscedasticity_test(
    wide,
    method = "nonparametric",
    imputation = "normal",
    seed = 77L
  )

  expect_equal(from_fit$statistic, from_wide$statistic, tolerance = 1e-10)
  expect_equal(from_fit$p.value, from_wide$p.value, tolerance = 1e-10)
  expect_equal(from_fit$pattern.counts, from_wide$pattern.counts)
})


test_that("mcar_homoscedasticity_test validates inputs", {
  expect_error(
    mcar_homoscedasticity_test(data.frame(a = 1:10, b = letters[1:10])),
    "must be numeric"
  )
  expect_error(
    mcar_homoscedasticity_test(matrix(1:20, ncol = 1)),
    "at least two variables"
  )
  expect_error(
    mcar_homoscedasticity_test(matrix(stats::rnorm(40), ncol = 2)),
    "requires missing values"
  )
  expect_error(
    mcar_homoscedasticity_test(make_jj_screening_data(), min_pattern_size = 1),
    "greater than or equal to 2"
  )
  expect_error(
    mcar_homoscedasticity_test(make_jj_screening_data(), nrep = 0),
    "positive integer"
  )
  expect_error(
    mcar_homoscedasticity_test(make_jj_screening_data(), n_min = 1),
    "greater than or equal to 2"
  )
  expect_error(
    mcar_homoscedasticity_test(make_jj_screening_data(), alpha = 1),
    "strictly between 0 and 1"
  )
  expect_error(
    mcar_homoscedasticity_test(make_jj_screening_data(), method = "invalid"),
    "'arg' should be one of"
  )
})
