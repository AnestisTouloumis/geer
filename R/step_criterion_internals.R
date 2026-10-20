compute_step_criterion <- function(object, criterion, cov_type) {
  values <- compute_criteria(
    object,
    cov_type = cov_type,
    criteria = criterion
  )
  as.numeric(values[[criterion]])
}


build_step_criterion_label <- function(criterion, cov_type) {
  if (criterion %in% c("QIC", "QICC")) {
    sprintf("%s with the %s covariance matrix", criterion, cov_type)
  } else {
    criterion
  }
}


compute_step_tolerance <- function(value) {
  sqrt(.Machine$double.eps) * max(1, abs(value))
}


build_step_criterion_key <- function(fitted_model, change) {
  build_step_model_key(
    stats::update.formula(
      stats::formula(fitted_model),
      stats::as.formula(paste(". ~ .", change))
    )
  )
}


fit_step_criterion_candidate <- function(fitted_model,
                                         change,
                                         criterion,
                                         cov_type,
                                         obs_no,
                                         cache) {
  key <- build_step_criterion_key(fitted_model, change)
  stored <- cache[[key]]
  if (!is.null(stored)) {
    return(stored)
  }
  fit <- update_step_fit(fitted_model, change, obs_no)
  candidate <- list(
    fit = fit,
    value = compute_step_criterion(fit, criterion, cov_type)
  )
  assign(key, candidate, envir = cache)
  candidate
}


build_step_candidate_table <- function(changes,
                                       candidates,
                                       fitted_model,
                                       current_value,
                                       criterion,
                                       step) {
  coef_no <- length(fitted_model$coefficients)
  aod <- data.frame(
    Df = c(
      NA_real_,
      vapply(
        candidates,
        function(x) abs(length(x$fit$coefficients) - coef_no),
        numeric(1)
      )
    ),
    value = c(current_value, vapply(candidates, `[[`, numeric(1), "value")),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(aod)[2L] <- criterion
  row.names(aod) <- c("<none>", changes)
  aod <- aod[order(aod[[criterion]], na.last = TRUE), , drop = FALSE]
  attr(aod, "heading") <- c(
    sprintf("Step %d candidates", step),
    "\nCurrent Model:",
    paste(
      deparse(stats::formula(fitted_model), width.cutoff = 500L),
      collapse = " "
    ),
    "\n"
  )
  class(aod) <- c("geer_anova", "anova", "data.frame")
  aod
}


evaluate_step_candidates <- function(fitted_model,
                                     fscope,
                                     direction,
                                     criterion,
                                     cov_type,
                                     obs_no,
                                     cache,
                                     current_value,
                                     step) {
  changes <- c(
    if (!identical(direction, "forward")) {
      vapply(
        fscope$drop,
        build_step_change_label,
        character(1),
        action = "drop"
      )
    },
    if (!identical(direction, "backward")) {
      vapply(
        fscope$add,
        build_step_change_label,
        character(1),
        action = "add"
      )
    }
  )
  if (!length(changes)) {
    return(NULL)
  }
  candidates <- lapply(
    changes,
    fit_step_criterion_candidate,
    fitted_model = fitted_model,
    criterion = criterion,
    cov_type = cov_type,
    obs_no = obs_no,
    cache = cache
  )
  values <- vapply(candidates, `[[`, numeric(1), "value")
  candidate_table <- build_step_candidate_table(
    changes = changes,
    candidates = candidates,
    fitted_model = fitted_model,
    current_value = current_value,
    criterion = criterion,
    step = step
  )
  ok <- is.finite(values)
  if (!any(ok)) {
    return(list(table = candidate_table, best = NULL))
  }
  best <- which(ok)[which.min(values[ok])]
  list(
    table = candidate_table,
    best = list(
      change = changes[[best]],
      fit = candidates[[best]]$fit,
      value = values[[best]]
    )
  )
}


build_step_criterion_results <- function(models,
                                         tables,
                                         fit,
                                         object,
                                         criterion,
                                         cov_type) {
  step_no <- vapply(models, `[[`, integer(1), "Step")
  aod <- data.frame(
    Step = ifelse(step_no == 0L, NA_integer_, step_no),
    Df = vapply(models, `[[`, NA_real_, "Df"),
    value = vapply(models, `[[`, NA_real_, "value"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(aod)[3L] <- criterion
  row.names(aod) <- build_step_row_labels(models)
  attr(aod, "heading") <- build_step_heading(
    object,
    fit,
    criterion_label = build_step_criterion_label(criterion, cov_type)
  )
  class(aod) <- c("geer_anova", "anova", "data.frame")
  fit$anova <- aod
  fit$candidates <- tables
  fit
}


run_step_criterion <- function(object,
                               scope,
                               direction,
                               criterion,
                               cov_type,
                               steps,
                               trace) {
  model_terms <- stats::terms(object)
  scope_factors <- extract_step_scope_factors(object, scope, model_terms)
  current_value <- compute_step_criterion(object, criterion, cov_type)
  if (!is.finite(current_value)) {
    stop(
      sprintf(
        "step_criterion failed: %s cannot be evaluated for the initial model",
        criterion
      ),
      call. = FALSE
    )
  }
  models <- list(list(
    Step = 0L,
    Change = NA_character_,
    Df = NA_real_,
    value = current_value
  ))
  tables <- list()
  cache <- new.env(parent = emptyenv())
  obs_no <- object$obs_no
  fitted_model <- object
  while (steps > 0L) {
    steps <- steps - 1L
    fscope <- stats::factor.scope(
      attr(model_terms, "factors"),
      list(add = scope_factors$add, drop = scope_factors$drop)
    )
    evaluation <- evaluate_step_candidates(
      fitted_model = fitted_model,
      fscope = fscope,
      direction = direction,
      criterion = criterion,
      cov_type = cov_type,
      obs_no = obs_no,
      cache = cache,
      current_value = current_value,
      step = length(models)
    )
    if (is.null(evaluation)) {
      break
    }
    tables[[length(tables) + 1L]] <- evaluation$table
    if (trace) {
      print(evaluation$table)
    }
    candidate <- evaluation$best
    threshold <- current_value - compute_step_tolerance(current_value)
    if (is.null(candidate) || !(candidate$value < threshold)) {
      break
    }
    df <- abs(
      length(candidate$fit$coefficients) - length(fitted_model$coefficients)
    )
    fitted_model <- candidate$fit
    model_terms <- stats::terms(fitted_model)
    current_value <- candidate$value
    models[[length(models) + 1L]] <- list(
      Step = length(models),
      Change = candidate$change,
      Df = df,
      value = current_value
    )
  }
  build_step_criterion_results(
    models = models,
    tables = tables,
    fit = fitted_model,
    object = object,
    criterion = criterion,
    cov_type = cov_type
  )
}
