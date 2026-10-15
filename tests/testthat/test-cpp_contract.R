cc_args <- function(...) {
  withr::local_seed(1)
  N <- 10
  a <- list(
    y_vector = rbinom(30, 1, 0.5), model_matrix = cbind(1, rnorm(30)),
    id_vector = rep(1:N, each = 3), repeated_vector = rep(1:3, N), weights_vector = rep(1, 30),
    link = "logit", family = "binomial", beta_vector = c(0, 0), offset = rep(0, 30),
    maxiter = 50L, tolerance = 1e-10, step_maxiter = 10L, step_multiplier = 1,
    jeffreys_power = 0.5, method = "gee", use_params = 1L, alpha_vector = numeric(0),
    alpha_fixed = 1L, correlation_structure = "independence", mdependence = 1L,
    phi = 1, phi_fixed = 1L, hold_nuisance = 0L
  )
  mod <- list(...)
  a[names(mod)] <- mod
  a
}

solve_cc <- function(...) do.call(fit_geesolver_cc, cc_args(...))

in_subprocess <- function(code) {
  dev <- requireNamespace("pkgload", quietly = TRUE) && isTRUE(pkgload::is_dev_package("geer"))
  if (dev) {
    out <- capture.output(eval(parse(text = code), envir = new.env(parent = asNamespace("geer"))))
    return(list(status = 0L, output = paste(out, collapse = "\n")))
  }
  old_libs <- Sys.getenv("R_LIBS", unset = NA)
  Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
  on.exit(if (is.na(old_libs)) Sys.unsetenv("R_LIBS") else Sys.setenv(R_LIBS = old_libs),
    add = TRUE
  )
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    c("-e", shQuote(paste0("suppressMessages(library(geer)); ", code))),
    stdout = TRUE, stderr = TRUE
  ))
  list(status = attr(out, "status") %||% 0L, output = paste(out, collapse = "\n"))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

or_solver_x <- c(
  -1.2, 0.3, 0.8, -0.5, 1.1, -0.2, 0.9, -1.4, 0.1, 0.4, -0.7, 1.3,
  -0.9, 0.6, -0.1, 1.5, -0.3, 0.2, -1.1, 0.7, -0.6, 0.0, 1.0, -0.8
)

or_solver_y <- c(
  0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1,
  0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0
)

call_or_solver <- function(y) {
  fit_geesolver_or(
    y, cbind(1, or_solver_x), rep(1:8, each = 3), rep(1:3, times = 8),
    rep(1, 24), "logit", c(0, 0), rep(0, 24), 100L, 1e-6, 10L, 1, 0.5,
    "gee", rep(1.5, 3)
  )
}

test_that("a missing association parameter is an R error, never a crash or an uninitialized read", {
  skip_on_cran()
  for (cs in c("exchangeable", "ar1")) {
    res <- in_subprocess(sprintf(
      "r <- tryCatch(geer:::get_correlation_matrix('%s', numeric(0), 3L), error = function(e) paste('R-ERROR:', conditionMessage(e))); cat(if (is.character(r)) r else 'NO-ERROR')", cs
    ))
    expect_identical(res$status, 0L, label = paste(cs, "process survived"))
    expect_match(res$output, "R-ERROR: .*(alpha|correlation)",
      label = paste(cs, "raised an R error naming the parameter")
    )
  }
  res <- in_subprocess("
    withr::local_seed(1); y <- rbinom(30, 1, .5); X <- cbind(1, rnorm(30))
    r <- tryCatch(geer:::fit_geesolver_cc(y, X, rep(1:10, each = 3), rep(1:3, 10), rep(1, 30), 'logit', 'binomial',
        c(0, 0), rep(0, 30), 50L, 1e-10, 10L, 1, .5, 'gee', 1L, numeric(0), 1L, 'exchangeable', 1L, 1, 1L, 0L),
        error = function(e) paste('R-ERROR:', conditionMessage(e))); cat(if (is.character(r)) r else 'NO-ERROR')")
  expect_identical(res$status, 0L)
  expect_match(res$output, "R-ERROR: .*(alpha|correlation)")
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
  expect_error(solve_cc(method = "nope"), "unsupported method")
  expect_error(solve_cc(correlation_structure = "nope"), "unsupported correlation structure")
  expect_error(solve_cc(link = "nope"))
  expect_error(solve_cc(family = "nope"))
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
