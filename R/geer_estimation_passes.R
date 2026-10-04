## Estimation passes shared by geewa(), geewa_binary() and the jackknife
## refits. The one-step methods are fitted in two passes:
##   bcgee-*         : GEE to convergence, then one brgee-* step;
##   hpgee-jeffreys  : independence PGEE to convergence, then one GEE step;
##   opgee-jeffreys  : independence PGEE to convergence, then one PGEE step.
## All other methods use a single pass.
##
## 'fit_pass(beta, pass, previous)' performs one call to the solver. 'pass' is
## a list with elements
##   method         : estimating method to use in this pass;
##   stage          : "single", "first" or "second";
##   independence   : TRUE when the pass uses an independence working structure;
##   one_step       : TRUE when only a single iteration is wanted;
##   carry_nuisance : TRUE when the association and dispersion parameters of
##                    the first pass are to be carried into this pass;
## and 'previous' is the first-pass fit (NULL unless stage is "second").
## 'check_first(fit, method)' must stop if the first pass is unusable.
run_geer_estimation_passes <- function(fit_pass, method, beta_start, check_first) {
  if (method %in% geer_bcgee_methods) {
    first_method <- "gee"
    second_method <- sub("bcgee", "brgee", method)
    independence <- FALSE
    carry_nuisance <- TRUE
  } else if (method %in% c("hpgee-jeffreys", "opgee-jeffreys")) {
    first_method <- "pgee-jeffreys"
    second_method <- if (identical(method, "hpgee-jeffreys")) "gee" else "pgee-jeffreys"
    independence <- TRUE
    carry_nuisance <- FALSE
  } else {
    pass <- list(method = method, stage = "single", independence = FALSE,
                 one_step = FALSE, carry_nuisance = FALSE)
    return(fit_pass(beta_start, pass, NULL))
  }
  first <- fit_pass(
    beta_start,
    list(method = first_method, stage = "first", independence = independence,
         one_step = FALSE, carry_nuisance = FALSE),
    NULL
  )
  check_first(first, method)
  fit_pass(
    as.numeric(first$beta_hat),
    list(method = second_method, stage = "second", independence = FALSE,
         one_step = TRUE, carry_nuisance = carry_nuisance),
    first
  )
}


## First-pass convergence check used by the main fitting functions.
check_geer_first_pass <- function(fit, method, tolerance) {
  last_criterion <- fit$criterion[ncol(fit$beta_mat) - 1L]
  if (last_criterion > tolerance) {
    if (method %in% geer_bcgee_methods) {
      stop("bias-corrected estimator is undefined because the corresponding GEE model did not converge", call. = FALSE)
    }
    stop(method, " estimator is undefined because the independence pgee-jeffreys model did not converge", call. = FALSE)
  }
  invisible(fit)
}
