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
    beta_start <- as.numeric(beta_start)
    p <- ncol(model_matrix)
    if (length(beta_start) != p) {
      stop("'beta_start' must be a numeric vector of length ", p, call. = FALSE)
    }
    return(beta_start)
  }
  control_glm <- do.call("brglm_control", control_glm)
  if (link != "identity" && !identical(family$family, "quasi")) {
    if (method %in% c("gee", "bcgee-naive", "bcgee-robust", "bcgee-empirical",
                      "brgee-robust", "brgee-empirical")) {
      type <- "ML"
    } else if (identical(method, "pgee-jeffreys")) {
      type <- "MPL_Jeffreys"
    } else {
      type <- "AS_mean"
    }

    glmfit <- try(
      brglmFit(
        x = model_matrix,
        y = y,
        family = family,
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
      ),
      silent = TRUE
    )
  } else {
    glmfit <- try(
      glm.fit(
        x = model_matrix,
        y = y,
        family = family,
        weights = weights,
        offset = offset,
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

  stop("cannot compute starting values; please supply 'beta_start'", call. = FALSE)
}

compute_geer_binary_start_values <- function(model_matrix,
                                             y,
                                             family,
                                             weights,
                                             offset,
                                             beta_start,
                                             control_glm,
                                             tolerance) {
  if (!is.null(beta_start)) {
    beta_start <- as.numeric(beta_start)
    p <- ncol(model_matrix)
    if (length(beta_start) != p) {
      stop("'beta_start' must be a numeric vector of length ", p, call. = FALSE)
    }
    return(beta_start)
  }
  control_glm <- do.call("brglm_control", control_glm)
  glmfit <- try(
    brglmFit(
      x = model_matrix,
      y = y,
      family = family,
      weights = weights,
      offset = offset,
      control = list(
        epsilon = tolerance,
        maxit = control_glm$maxit,
        type = "AS_mean",
        trace = FALSE,
        slowit = control_glm$slowit,
        max_step_factor = control_glm$max_step_factor
      )
    ),
    silent = TRUE
  )
  if (!inherits(glmfit, "try-error") && all(is.finite(glmfit$coefficients))) {
    return(glmfit$coefficients)
  }
  stop("cannot compute starting values; please supply 'beta_start'", call. = FALSE)
}
