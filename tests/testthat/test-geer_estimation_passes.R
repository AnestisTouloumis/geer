test_that("run_geer_estimation_passes reproduces the pass schedule per method", {
  run <- function(method) {
    calls <- list()
    fit_pass <- function(beta, pass, previous) {
      calls[[length(calls) + 1L]] <<- list(beta = beta, pass = pass,
                                           has_previous = !is.null(previous))
      list(beta_hat = beta + 1, alpha = 0.5, phi = 2)
    }
    checked <- character()
    out <- run_geer_estimation_passes(
      fit_pass, method, beta_start = c(0, 0),
      check_first = function(fit, method) checked <<- c(checked, method)
    )
    list(calls = calls, checked = checked, out = out)
  }

  single <- run("brgee-robust")
  expect_length(single$calls, 1L)
  expect_identical(single$calls[[1L]]$pass$method, "brgee-robust")
  expect_identical(single$calls[[1L]]$pass$stage, "single")
  expect_length(single$checked, 0L)

  bc <- run("bcgee-robust")
  expect_equal(vapply(bc$calls, function(z) z$pass$method, ""),
               c("gee", "brgee-robust"))
  expect_true(bc$calls[[2L]]$pass$carry_nuisance)
  expect_true(bc$calls[[2L]]$pass$one_step)
  expect_false(bc$calls[[1L]]$pass$independence)
  expect_equal(bc$calls[[2L]]$beta, c(1, 1))
  expect_identical(bc$checked, "bcgee-robust")

  hp <- run("hpgee-jeffreys")
  expect_equal(vapply(hp$calls, function(z) z$pass$method, ""),
               c("pgee-jeffreys", "gee"))
  expect_true(hp$calls[[1L]]$pass$independence)
  expect_false(hp$calls[[2L]]$pass$independence)
  expect_false(hp$calls[[2L]]$pass$carry_nuisance)

  op <- run("opgee-jeffreys")
  expect_equal(vapply(op$calls, function(z) z$pass$method, ""),
               c("pgee-jeffreys", "pgee-jeffreys"))
})

test_that("check_geer_first_pass stops only on non-convergence", {
  fit <- list(criterion = c(1, 1e-10, 0), beta_mat = matrix(0, 2, 3))
  expect_silent(check_geer_first_pass(fit, "bcgee-naive", 1e-6))
  bad <- list(criterion = c(1, 0.5, 0), beta_mat = matrix(0, 2, 3))
  expect_error(check_geer_first_pass(bad, "bcgee-naive", 1e-6),
               "bias-corrected estimator is undefined")
  expect_error(check_geer_first_pass(bad, "hpgee-jeffreys", 1e-6),
               "hpgee-jeffreys estimator is undefined")
})
