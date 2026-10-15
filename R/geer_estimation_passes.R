run_geer_estimation_passes <- function(fit_pass,
                                       method,
                                       beta_start,
                                       check_first) {
  if (method %in% geer_bcgee_methods) {
    first_method <- "gee"
    second_method <- sub("bcgee", "brgee", method)
    independence <- FALSE
    carry_nuisance <- TRUE
    hold_nuisance <- FALSE
  } else if (method %in% c("hpgee-jeffreys", "opgee-jeffreys")) {
    first_method <- "pgee-jeffreys"
    second_method <- if (identical(method, "hpgee-jeffreys")) {
      "gee"
    } else {
      "pgee-jeffreys"
    }
    independence <- TRUE
    carry_nuisance <- FALSE
    hold_nuisance <- TRUE
  } else {
    pass <- list(
      method = method,
      stage = "single",
      independence = FALSE,
      one_step = FALSE,
      carry_nuisance = FALSE,
      hold_nuisance = FALSE
    )
    return(fit_pass(beta_start, pass, NULL))
  }
  first <- fit_pass(
    beta_start,
    list(
      method = first_method,
      stage = "first",
      independence = independence,
      one_step = FALSE,
      carry_nuisance = FALSE,
      hold_nuisance = FALSE
    ),
    NULL
  )
  check_first(first, method)
  fit_pass(
    as.numeric(first$beta_hat),
    list(
      method = second_method,
      stage = "second",
      independence = FALSE,
      one_step = TRUE,
      carry_nuisance = carry_nuisance,
      hold_nuisance = hold_nuisance
    ),
    first
  )
}


get_solver_failure <- function(fit) {
  failure <- fit$failure
  if (is.character(failure) && length(failure) == 1L && !is.na(failure)) {
    failure
  } else {
    ""
  }
}


check_geer_first_pass <- function(fit, method, tolerance) {
  last_criterion <- fit$criterion[ncol(fit$beta_mat) - 1L]
  if (!isTRUE(last_criterion <= tolerance)) {
    reason <- get_solver_failure(fit)
    detail <- if (nzchar(reason)) paste0(": ", reason) else ""
    if (method %in% geer_bcgee_methods) {
      stop(
        "bias-corrected estimator is undefined because the corresponding GEE ",
        "model did not converge",
        detail,
        call. = FALSE
      )
    }
    stop(
      method,
      " estimator is undefined because the independence pgee-jeffreys model ",
      "did not converge",
      detail,
      call. = FALSE
    )
  }
  invisible(fit)
}
