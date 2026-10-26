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


check_geer_first_pass <- function(fit, method, tolerance, fit_function) {
  last_criterion <- fit$criterion[ncol(fit$beta_mat) - 1L]
  if (!isTRUE(last_criterion <= tolerance)) {
    reason <- get_solver_failure(fit)
    detail <- if (nzchar(reason)) paste0(": ", reason) else ""
    first_stage <- if (method %in% geer_bcgee_methods) {
      "GEE"
    } else {
      "independence pgee-jeffreys"
    }
    stop(
      sprintf(
        paste0(
          "%s failed: the %s estimator is undefined because its first-stage ",
          "%s fit did not converge%s"
        ),
        fit_function,
        method,
        first_stage,
        detail
      ),
      call. = FALSE
    )
  }
  invisible(fit)
}


run_geer_solver <- function(expr, context) {
  check_choice(context, c("geewa", "geewa_binary"), "context")
  tryCatch(
    force(expr),
    error = function(error) {
      message <- conditionMessage(error)
      if (!inherits(error, "GeerInputError")) {
        message <- sub("^[^:]+ failed: ", "", message)
        message <- paste0(context, " failed: ", message)
      }
      stop(message, call. = FALSE)
    }
  )
}
