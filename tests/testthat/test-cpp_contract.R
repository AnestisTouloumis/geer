list2env(get_cached_test_fixtures("cpp_contract"), environment())

test_that("a missing association parameter is an R error, never a crash or an uninitialized read", {
  skip_on_cran()
  res <- in_subprocess("
    report <- function(label, expr) {
      r <- tryCatch(
        expr,
        error = function(e) paste('R-ERROR:', conditionMessage(e))
      )
      cat(label, if (is.character(r)) r else 'NO-ERROR', '| ')
    }
    gc <- geer:::get_correlation_matrix_cc
    report('exchangeable', gc('exchangeable', numeric(0), 3L))
    report('ar1', gc('ar1', numeric(0), 3L))
    withr::local_seed(1)
    y <- rbinom(30, 1, .5)
    X <- cbind(1, rnorm(30))
    report('solver', geer:::fit_geesolver_cc(
      y, X, rep(1:10, each = 3), rep(1:3, 10), rep(1, 30), 'logit',
      'binomial', c(0, 0), rep(0, 30), 50L, 1e-10, 10L, 1, .5, 'gee', 1L,
      numeric(0), 1L, 'exchangeable', 1L, 1, 1L, 0L
    ))")
  expect_identical(res$status, 0L)
  for (label in c("exchangeable", "ar1", "solver")) {
    expect_match(
      res$output,
      paste0(label, " R-ERROR: [^|]*(alpha|correlation)"),
      label = label
    )
  }
})

test_that("observations of a cluster must be contiguous: interleaved ids are rejected", {
  ok <- solve_cc()
  o <- order(rep(1:3, 10) + 0.001 * seq_len(30))
  a <- cc_args()
  bad <- function() {
    do.call(
      fit_geesolver_cc,
      modifyList(a, list(
        y_vector = a$y_vector[o], model_matrix = a$model_matrix[o, ],
        id_vector = a$id_vector[o], repeated_vector = a$repeated_vector[o]
      ))
    )
  }
  expect_error(bad())
})

test_that("inconsistent argument lengths and empty inputs are errors", {
  expect_error(solve_cc(weights_vector = rep(1, 29)))
  expect_error(solve_cc(beta_vector = c(0, 0, 0)))
  expect_error(solve_cc(offset = rep(0, 29)))
  expect_error(solve_cc(
    y_vector = numeric(0), model_matrix = matrix(0, 0, 2), id_vector = numeric(0),
    repeated_vector = numeric(0), weights_vector = numeric(0), offset = numeric(0)
  ))
  expect_error(get_estimating_equations_cc(
    numeric(0), matrix(0, 0, 2), numeric(0), numeric(0), numeric(0),
    "logit", "binomial", numeric(0), numeric(0), "independence",
    numeric(0), 1
  ), "empty")
})

test_that("invalid time indices are errors, not silent mis-indexing", {
  for (bad in list(NA_real_, 0, -1, 1.5)) {
    rep_ <- rep(1:3, 10)
    rep_[2] <- bad
    expect_error(
      solve_cc(
        repeated_vector = rep_, correlation_structure = "exchangeable",
        alpha_vector = 0.2
      ),
      label = paste("repeated =", bad)
    )
  }
})

test_that("non-finite response, offset and weights are errors at the C++ boundary", {
  y <- cc_args()$y_vector
  y[3] <- NaN
  expect_error(solve_cc(y_vector = y))
  expect_error(solve_cc(offset = replace(rep(0, 30), 3, NaN)))
  expect_error(solve_cc(weights_vector = replace(rep(1, 30), 3, NaN)))
  expect_error(solve_cc(phi = 0))
  expect_error(solve_cc(phi = NaN))
})

test_that("unknown method / structure / link / family names are errors", {
  expect_error(solve_cc(method = "nope"), "'method' must be one of")
  expect_error(
    solve_cc(correlation_structure = "nope"),
    "'correlation_structure' must be one of"
  )
  expect_error(solve_cc(link = "nope"), "'link' must be one of")
  expect_error(solve_cc(family = "nope"), "'family' must be one of")
})

test_that("negative iteration limits are rejected with an informative error", {
  expect_error(solve_cc(maxiter = -1L), "maxiter")
  expect_error(solve_cc(step_maxiter = -1L), "maxiter")
})

test_that("the odds-ratio solver fits binary and proportion responses", {
  binary <- call_or_solver(or_solver_y)
  expect_identical(binary$failure, "")
  expect_true(all(is.finite(binary$beta_hat)))
  iterations <- sum(binary$criterion != 0)
  expect_lte(binary$criterion[iterations], 1e-6)

  proportions <- or_solver_y
  proportions[c(3, 8)] <- c(0.3, 0.75)
  fractional <- call_or_solver(proportions)
  expect_identical(fractional$failure, "")
  expect_true(all(is.finite(fractional$beta_hat)))
})

test_that("the odds-ratio solver rejects responses outside [0, 1]", {
  bad_above <- or_solver_y
  bad_above[4] <- 2
  expect_error(
    call_or_solver(bad_above),
    "the response must lie in \\[0, 1\\]; observation 4 has value 2"
  )

  bad_below <- or_solver_y
  bad_below[11] <- -1
  expect_error(call_or_solver(bad_below), "observation 11 has value -1")

  for (bad in c(NaN, Inf, -Inf)) {
    bad_value <- or_solver_y
    bad_value[6] <- bad
    expect_error(call_or_solver(bad_value), "observation 6 has value")
  }
})

test_that("indefinite working covariances stop at both boundaries", {
  expect_no_warning(
    expect_error(
      get_covariance_matrices_or(
        c(0, 1, 0), matrix(c(-1, 1, 1), 3L, 1L), rep(1, 3L), 1:3,
        rep(1, 3L), "logit", rep(0.5, 3L), rep(0, 3L), c(50, 50, 0.02)
      ),
      paste0(
        "get_covariance_matrices_or failed: cluster 1 (id = 1): ",
        "get_v_matrix_or failed: the working covariance matrix is not ",
        "positive definite"
      ),
      fixed = TRUE
    )
  )
  expect_no_warning(
    expect_error(
      geewa_binary(
        formula = y ~ 1,
        data = data.frame(y = c(0, 1, 0), id = 1, visit = 1:3),
        id = id, repeated = visit, orstr = "fixed", beta_start = 0,
        alpha_vector = c(50, 50, 0.02)
      ),
      paste0(
        "^geewa_binary failed:.*working covariance matrix is not ",
        "positive definite"
      )
    )
  )
})


test_that("the odds-ratio boundary validates iterations and finite inputs", {
  args <- list(
    y_vector = or_solver_y, model_matrix = cbind(1, or_solver_x),
    id_vector = rep(1:8, each = 3L), repeated_vector = rep(1:3, 8L),
    weights_vector = rep(1, 24L), link = "logit", beta_vector = c(0, 0),
    offset = rep(0, 24L), maxiter = 100L, tolerance = 1e-6,
    step_maxiter = 10L, step_multiplier = 1, jeffreys_power = 0.5,
    method = "gee", alpha_vector = rep(1.5, 3L)
  )
  for (field in c("maxiter", "step_maxiter")) {
    invalid <- args
    invalid[[field]] <- -1L
    expect_error(do.call(fit_geesolver_or, invalid))
  }
  for (field in c("model_matrix", "offset", "weights_vector", "beta_vector")) {
    invalid <- args
    invalid[[field]][1L] <- NaN
    expect_error(do.call(fit_geesolver_or, invalid))
  }
})
