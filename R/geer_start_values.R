build_family_no_aic <- function(family) {
  family$aic <- function(y, n, mu, wt, dev) NA_real_
  family
}


evaluate_noninteger_quiet <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("^non-integer (x|#successes)", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}


check_beta_start <- function(beta_start, p) {
  if (!is.numeric(beta_start) || length(beta_start) != p) {
    stop("'beta_start' must be a numeric vector of length ", p, call. = FALSE)
  }
  if (any(!is.finite(beta_start))) {
    stop("'beta_start' must contain only finite values", call. = FALSE)
  }
  as.numeric(beta_start)
}


normalize_control_glm <- function(control_glm) {
  if (!is.list(control_glm)) {
    stop("'control_glm' must be a list", call. = FALSE)
  }
  do.call("brglm_control", control_glm)
}


get_start_type <- function(method) {
  if (method %in% c("gee", geer_bcgee_methods)) {
    return("ML")
  }
  if (identical(method, "brgee-naive")) {
    return("AS_mean")
  }
  "MPL_Jeffreys"
}


compute_geer_start_values <- function(model_matrix,
                                      y,
                                      family,
                                      weights,
                                      offset,
                                      method,
                                      link,
                                      beta_start,
                                      control,
                                      control_glm) {
  if (!is.null(beta_start)) {
    return(check_beta_start(beta_start, ncol(model_matrix)))
  }
  control_glm <- normalize_control_glm(control_glm)
  if (link != "identity" &&
      !(family$family %in% c("quasi", "quasibinomial", "quasipoisson"))) {
    type <- get_start_type(method)
    glmfit <- try(
      evaluate_noninteger_quiet(brglmFit(
        x = model_matrix,
        y = y,
        family = build_family_no_aic(family),
        weights = weights,
        offset = offset,
        control = list(
          epsilon = control$tolerance,
          maxit = control_glm$maxit,
          type = type,
          trace = FALSE,
          slowit = control_glm$slowit,
          max_step_factor = control_glm$max_step_factor,
          a = control$jeffreys_power
        )
      )),
      silent = TRUE
    )
  } else {
    glmfit <- try(
      stats::glm.fit(
        x = model_matrix,
        y = y,
        weights = weights,
        offset = offset,
        family = family,
        control = list(
          epsilon = control$tolerance,
          maxit = control_glm$maxit,
          trace = FALSE
        )
      ),
      silent = TRUE
    )
  }
  if (!inherits(glmfit, "try-error") && all(is.finite(glmfit$coefficients))) {
    return(glmfit$coefficients)
  }
  stop(
    "starting values could not be computed; supply 'beta_start'",
    call. = FALSE
  )
}


compute_binary_start_values <- function(model_matrix,
                                        y,
                                        family,
                                        weights,
                                        offset,
                                        method,
                                        beta_start,
                                        control_glm,
                                        tolerance,
                                        jeffreys_power) {
  if (!is.null(beta_start)) {
    return(check_beta_start(beta_start, ncol(model_matrix)))
  }
  type <- get_start_type(method)
  control_glm <- normalize_control_glm(control_glm)
  glmfit <- try(
    evaluate_noninteger_quiet(brglmFit(
      x = model_matrix,
      y = y,
      family = build_family_no_aic(family),
      weights = weights,
      offset = offset,
      control = list(
        epsilon = tolerance,
        maxit = control_glm$maxit,
        type = type,
        trace = FALSE,
        slowit = control_glm$slowit,
        max_step_factor = control_glm$max_step_factor,
        a = jeffreys_power
      )
    )),
    silent = TRUE
  )
  if (!inherits(glmfit, "try-error") && all(is.finite(glmfit$coefficients))) {
    return(glmfit$coefficients)
  }
  stop(
    "starting values could not be computed; supply 'beta_start'",
    call. = FALSE
  )
}
