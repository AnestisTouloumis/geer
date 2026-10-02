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
    neyman_nulls = geer:::jj_neyman_nulls(c(8L, 8L, 8L), nrep = 20L, n_min = 2L)
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
    "only 8 complete cases are available among the retained rows"
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
    expect_warning(
      out <- mcar_homoscedasticity_test(
        x,
        method = "nonparametric",
        imputation = "normal",
        seed = 9L
      ),
      "can inflate the size of the nonparametric test"
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

  expect_warning(
    from_fit <- mcar_homoscedasticity_test(
      fit,
      method = "nonparametric",
      imputation = "normal",
      seed = 77L
    ),
    "can inflate the size of the nonparametric test"
  )
  expect_warning(
    from_wide <- mcar_homoscedasticity_test(
      wide,
      method = "nonparametric",
      imputation = "normal",
      seed = 77L
    ),
    "can inflate the size of the nonparametric test"
  )

  expect_equal(from_fit$statistic, from_wide$statistic, tolerance = 1e-10)
  expect_equal(from_fit$p.value, from_wide$p.value, tolerance = 1e-10)
  expect_equal(from_fit$pattern.counts, from_wide$pattern.counts)
  expect_identical(from_fit$data.name, "repeated response from fitted geer object")
  expect_identical(from_wide$data.name, "wide")
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
  expect_error(
    mcar_homoscedasticity_test(make_jj_screening_data(), n_imputations = 0),
    "'n_imputations' must be a positive integer"
  )
  expect_error(
    mcar_homoscedasticity_test(make_jj_screening_data(), seed = 1e10),
    "'seed' must be NULL or a single whole number"
  )
  expect_error(
    mcar_homoscedasticity_test(make_jj_screening_data(), seed = 2.5),
    "'seed' must be NULL or a single whole number"
  )
  expect_error(
    mcar_homoscedasticity_test(make_jj_screening_data(), nrep = "200"),
    "'nrep' must be a positive integer"
  )
  expect_error(
    mcar_homoscedasticity_test(make_jj_screening_data(), nrep = TRUE),
    "'nrep' must be a positive integer"
  )
  expect_error(
    mcar_homoscedasticity_test(make_jj_screening_data(), n_min = "2"),
    "'n_min' must be an integer greater than or equal to 2"
  )
  expect_error(
    mcar_homoscedasticity_test(make_jj_screening_data(), maxit = 1e12),
    "'maxit' must be a positive integer"
  )
  expect_error(
    mcar_homoscedasticity_test(make_jj_screening_data(), alpha = "0.05"),
    "'alpha' must be a single finite numeric value"
  )
  expect_error(
    mcar_homoscedasticity_test(make_jj_screening_data(), tol = "1e-8"),
    "'tol' must be a single positive finite number"
  )
  expect_error(
    mcar_homoscedasticity_test(make_jj_screening_data(), data = data.frame()),
    "'data' is only used when 'object' is a fitted 'geer' object"
  )
})


test_that("Hawkins transformed statistics match MissMech on deterministic data", {
  i <- 1:30
  completed <- cbind(
    sin(i),
    cos(1.7 * i) + 0.3 * sin(i),
    sin(0.5 * i) * cos(i)
  )
  group <- rep(1:3, c(12L, 10L, 8L))
  counts <- c(12L, 10L, 8L)

  # Reference values from MissMech 1.0.4: Hawkins(), TestUNey()$n4 and
  # AndersonDarling() applied to the same completed data.
  missmech_fij <- c(
    0.5491237911, 1.2542372814, 1.4383803944, 1.3198396708, 0.9390506887,
    0.4078911978, 0.8520623071, 0.9178865958, 1.8588686954, 1.1126646308,
    1.2270582680, 0.2333683770, 0.7087359739, 0.4812351325, 1.5367082060,
    1.5844134895, 1.6246042483, 0.9496881110, 0.2639144208, 0.7883786779,
    0.6697750429, 2.2476750710, 0.8148038203, 1.2453297129, 0.0167763106,
    1.7747780658, 0.9693748385, 1.6217738517, 0.6906485960, 0.7873700661
  )

  hawkins <- geer:::jj_hawkins_test(
    completed,
    group = group,
    group_counts = counts,
    neyman_nulls = geer:::jj_neyman_nulls(counts, nrep = 10L, n_min = 2L)
  )
  expect_equal(hawkins$f.values, missmech_fij, tolerance = 1e-9)
  expect_equal(hawkins$group.statistics[1L], 5.7707840618, tolerance = 1e-9)

  ad <- geer:::jj_anderson_darling_test(hawkins$f.values, group, counts)
  expect_equal(ad$statistic, 0.8206642148, tolerance = 1e-9)
  expect_equal(ad$variance, 0.9910839542, tolerance = 1e-9)
  expect_equal(
    ad$group.statistics,
    c(0.2930375607, 0.3227273003, 0.2048993537),
    tolerance = 1e-9
  )
  expect_equal(sum(ad$group.statistics), ad$statistic, tolerance = 1e-12)
})


test_that("Anderson-Darling statistic reproduces the published MissMech example", {
  set.seed(50)
  x <- c(stats::rnorm(30), stats::runif(45), stats::rnorm(60, 2, 3))
  counts <- c(30L, 45L, 60L)
  out <- geer:::jj_anderson_darling_test(x, rep(1:3, counts), counts)

  # Jamshidian, Jalal and Jansen (2014, Section 5.1).
  expect_equal(out$statistic, 18.62008, tolerance = 1e-6)
  expect_equal(out$variance, 1.119502, tolerance = 1e-6)
  expect_equal(
    out$group.statistics,
    c(6.566425, 5.075349, 6.978307),
    tolerance = 1e-6
  )
  expect_true(out$extrapolated)
})


test_that("Anderson-Darling p-values match kSamples reference quantiles", {
  # Reference values from kSamples 1.2-12: ad.pval(standardized, m, version = 1).
  reference <- data.frame(
    m = rep(c(1, 2, 4, 9), times = 5L),
    standardized = rep(c(-1.5, 0.5, 2, 3.5, 6), each = 4L),
    p.value = c(
      1, 0.999800276218, 0.986828652144, 0.965962313201,
      0.20763624818, 0.234448984674, 0.255492542308, 0.272532810456,
      0.0483861313739, 0.0470203946581, 0.0439580739493, 0.0390517529971,
      0.0127506350622, 0.009382004813, 0.00627451403749, 0.00366513941818,
      0.00151990596119, 0.000617335256168, 0.00019898879606, 4.76351735552e-05
    )
  )

  p_values <- mapply(
    function(standardized, m) geer:::jj_ad_p_value(standardized, m)$p.value,
    reference$standardized,
    reference$m
  )
  expect_equal(p_values, reference$p.value, tolerance = 1e-6)

  expect_identical(dim(geer:::jj_ad_quantiles), c(35L, 8L))
  expect_equal(
    geer:::jj_ad_quantiles[1L, ],
    c(-1.1954, -1.5806, -1.8172, -2.0032, -2.2526, -2.4204, -2.5283, -4.2649)
  )
  expect_equal(
    geer:::jj_ad_quantiles[35L, ],
    c(11.8537, 9.5482, 8.5568, 8.0283, 7.4418, 6.9524, 6.6195, 4.2649)
  )
  expect_false(geer:::jj_ad_p_value(2, 4)$extrapolated)
})


test_that("simulated Neyman p-values use the (1 + b) / (nrep + 1) form", {
  set.seed(31)
  x <- stats::runif(12)
  null <- geer:::jj_neyman_null(12L, nrep = 499L)
  out <- geer:::jj_neyman_p_value(x, null = null)

  expect_true(out$simulated)
  expect_equal(
    out$p.value,
    (1 + sum(null >= geer:::jj_neyman_statistic(x))) / 500
  )
  expect_gte(out$p.value, 1 / 500)

  set.seed(8)
  looped <- vapply(
    seq_len(25L),
    function(i) geer:::jj_neyman_statistic(stats::runif(7L)),
    numeric(1)
  )
  set.seed(8)
  blocked <- geer:::jj_neyman_null(7L, nrep = 25L, block_size = 10L)
  expect_equal(blocked, looped, tolerance = 1e-12)
})


test_that("Hawkins transformation stops when a case deletion is singular", {
  completed <- rbind(
    c(0, 0), c(1, 1), c(2, 2),
    c(0, 0), c(1, 1), c(2, 2.5),
    c(0, 1), c(1, 2), c(2, 3)
  )
  expect_error(
    geer:::jj_hawkins_test(
      completed,
      group = rep(1:3, each = 3L),
      group_counts = c(3L, 3L, 3L),
      neyman_nulls = NULL,
      test_uniformity = FALSE
    ),
    "leaves a singular pooled covariance matrix"
  )
})


test_that("multiple imputations are analyzed and summarized", {
  x <- make_jj_screening_data()
  single <- mcar_homoscedasticity_test(
    x,
    method = "auto",
    n_min = 2L,
    seed = 44L
  )
  multiple <- mcar_homoscedasticity_test(
    x,
    method = "auto",
    n_imputations = 4L,
    n_min = 2L,
    seed = 44L
  )

  expect_identical(multiple$n.imputations, 4L)
  expect_identical(nrow(multiple$imputations$tests), 4L)
  expect_identical(dim(multiple$imputations$hawkins.group.p.values), c(4L, 3L))
  expect_identical(
    dim(multiple$imputations$nonparametric.group.statistics),
    c(4L, 3L)
  )
  expect_identical(
    colnames(multiple$imputations$hawkins.group.p.values),
    names(multiple$pattern.counts)
  )
  expect_equal(multiple$imputed.data, single$imputed.data)
  expect_equal(multiple$p.value, single$p.value)
  expect_equal(
    multiple$imputations$tests$hawkins.p.value[1L],
    multiple$hawkins$p.value
  )
  expect_equal(
    multiple$imputations$tests$nonparametric.p.value[1L],
    multiple$nonparametric$p.value
  )
  expect_true(length(unique(multiple$imputations$tests$nonparametric.statistic)) > 1L)
})


test_that("user-supplied completed data bypass imputation", {
  set.seed(1010)
  y1 <- matrix(stats::rnorm(50 * 3), ncol = 3)
  y2 <- matrix(stats::rnorm(40 * 3), ncol = 3)
  y3 <- matrix(stats::rnorm(30 * 3, sd = sqrt(2)), ncol = 3)
  complete <- rbind(y1, y2, y3)
  incomplete <- complete
  incomplete[1:50, 1] <- NA_real_
  incomplete[51:90, 2] <- NA_real_
  incomplete[91:120, 3] <- NA_real_

  out <- mcar_homoscedasticity_test(
    incomplete,
    method = "nonparametric",
    imputed_data = complete
  )
  direct <- geer:::jj_anderson_darling_test(
    geer:::jj_hawkins_test(
      complete,
      group = rep(1:3, c(50L, 40L, 30L)),
      group_counts = c(50L, 40L, 30L),
      neyman_nulls = NULL,
      test_uniformity = FALSE
    )$f.values,
    rep(1:3, c(50L, 40L, 30L)),
    c(50L, 40L, 30L)
  )

  expect_identical(out$imputation.used, "user-supplied")
  expect_null(out$location)
  expect_equal(unname(out$imputed.data), complete)
  expect_equal(unname(out$statistic), direct$statistic)
  expect_equal(out$p.value, direct$p.value)

  expect_error(
    mcar_homoscedasticity_test(incomplete, imputed_data = complete[-1L, ]),
    "same dimensions"
  )
  perturbed <- complete
  perturbed[60L, 1L] <- perturbed[60L, 1L] + 1
  expect_error(
    mcar_homoscedasticity_test(incomplete, imputed_data = perturbed),
    "reproduce the observed values"
  )
  with_na <- complete
  with_na[1L, 1L] <- NA_real_
  expect_error(
    mcar_homoscedasticity_test(incomplete, imputed_data = with_na),
    "missing or non-finite"
  )
  expect_error(
    mcar_homoscedasticity_test(
      incomplete,
      n_imputations = 2L,
      imputed_data = complete
    ),
    "must be 1 when 'imputed_data' is supplied"
  )
  expect_error(
    mcar_homoscedasticity_test(
      incomplete,
      imputation = "normal",
      imputed_data = complete
    ),
    "'imputation' must not be supplied when 'imputed_data' is supplied"
  )
  expect_error(
    mcar_homoscedasticity_test(
      incomplete,
      imputation = "distribution-free",
      imputed_data = complete
    ),
    "'imputation' must not be supplied when 'imputed_data' is supplied"
  )
})

test_that("jj_hawkins_test requires one Neyman null entry per pattern", {
  x <- make_jj_screening_data()
  x[is.na(x)] <- 0
  expect_error(
    geer:::jj_hawkins_test(
      x,
      group = rep(1:3, c(30L, 15L, 15L)),
      group_counts = c(30L, 15L, 15L),
      neyman_nulls = list(NULL, NULL)
    ),
    "'neyman_nulls' must be a list with one element per retained missingness pattern"
  )
})
