.step_p_scope_factors <- function(object, scope, model_terms) {
  if (!missing(scope) && !is.null(scope) &&
      !is.list(scope) && !inherits(scope, "formula") && !is.character(scope)) {
    stop("'scope' must be NULL, a formula, a character vector, or a list", call. = FALSE)
  }
  if (missing(scope)) {
    factors_drop <- numeric()
    factors_add <- attr(model_terms, "factors")
  } else if (is.list(scope)) {
    factors_drop <- scope$lower
    factors_add <- scope$upper

    factors_drop <- if (!is.null(factors_drop)) {
      attr(stats::terms(stats::update.formula(object, factors_drop)), "factors")
    } else {
      numeric()
    }

    factors_add <- if (!is.null(factors_add)) {
      attr(stats::terms(stats::update.formula(object, factors_add)), "factors")
    } else {
      numeric()
    }
  } else {
    factors_add <- if (!is.null(scope)) {
      attr(stats::terms(stats::update.formula(object, scope)), "factors")
    } else {
      numeric()
    }
    factors_drop <- numeric()
  }

  list(
    add = factors_add,
    drop = factors_drop
  )
}


.step_p_init_models <- function(object, cov_type, steps) {
  steps <- as.integer(steps)

  models <- vector("list", steps + 1L)
  models[[1L]] <- list(
    Step = "",
    Df = NA_real_,
    Chi = NA_real_,
    `Pr(>Chi)` = NA_real_,
    CIC = extract_cic(object, cov_type)
  )

  list(
    steps = steps,
    models = models,
    models_no = 1L
  )
}


.step_p_update_fit <- function(fitted_model, change, obs_no) {
  fitted_model <- stats::update(
    fitted_model,
    formula = paste(". ~ .", change),
    evaluate = FALSE
  )
  fitted_model <- eval.parent(fitted_model)

  obs_no_new <- fitted_model$obs_no
  if (all(is.finite(c(obs_no, obs_no_new))) && obs_no_new != obs_no) {
    stop("number of rows in use has changed: remove missing values?", call. = FALSE)
  }

  fitted_model
}


.step_p_append_model <- function(models, models_no, step, df, chi, pval, cic) {
  models_no <- models_no + 1L

  models[[models_no]] <- list(
    Step = step,
    Df = df,
    Chi = chi,
    `Pr(>Chi)` = pval,
    CIC = cic
  )

  list(
    models = models,
    models_no = models_no
  )
}



## printing results
.step_p_results <- function(models, fit, object) {
  aod <- data.frame(
    Step = vapply(models, `[[`, "", "Step"),
    Df = vapply(models, `[[`, NA_real_, "Df"),
    Chi = vapply(models, `[[`, NA_real_, "Chi"),
    `Pr(>Chi)` = vapply(models, `[[`, NA_real_, "Pr(>Chi)"),
    CIC = vapply(models, `[[`, NA_real_, "CIC"),
    check.names = FALSE
  )
  heading <- c(
    "Stepwise Model Path \nAnalysis of Deviance Table",
    "\nInitial Model:", deparse(stats::formula(object)),
    "\nFinal Model:", deparse(stats::formula(fit)), "\n"
  )
  attr(aod, "heading") <- heading
  class(aod) <- c("anova", "data.frame")
  fit$anova <- aod
  fit
}

.step_p_selected_row_backward <- function(aod) {
  aod2 <- aod[-1L, , drop = FALSE]
  if (!nrow(aod2)) {
    return(NA_integer_)
  }

  zdf <- !is.na(aod2$Df) & aod2$Df == 0
  idx <- which(zdf)
  if (length(idx)) {
    return(idx[length(idx)] + 1L)
  }


  pvals <- aod2[, "Pr(>Chi)"]
  ok <- is.finite(pvals)
  if (!any(ok)) {
    return(NA_integer_)
  }

  which(ok)[which.max(pvals[ok])] + 1L
}



.step_p_selected_row_forward <- function(aod) {
  aod2 <- aod[-1L, , drop = FALSE]
  if (!nrow(aod2)) {
    return(NA_integer_)
  }

  pvals <- aod2[, "Pr(>Chi)"]
  ok <- is.finite(pvals)
  if (!any(ok)) {
    return(NA_integer_)
  }

  which(ok)[which.min(pvals[ok])] + 1L
}

.step_p_backward_should_stop <- function(aod, row, pvalue) {
  if (is.na(row) || row <= 1L || row > nrow(aod)) {
    return(TRUE)
  }

  aod2 <- aod[-1L, , drop = FALSE]
  zdf <- !is.na(aod2$Df) & aod2$Df == 0
  if (any(zdf)) {
    return(FALSE)
  }

  pval <- aod[row, "Pr(>Chi)"]
  is.na(pval) || pval < pvalue
}


.step_p_forward_should_stop <- function(aod, row, pvalue) {
  if (is.na(row) || row <= 1L || row > nrow(aod)) {
    return(TRUE)
  }

  pval <- aod[row, "Pr(>Chi)"]
  is.na(pval) || pval > pvalue
}

.step_p_selected_step_label <- function(aod, row) {
  if (is.na(row) || row < 1L || row > nrow(aod)) {
    return(NA_character_)
  }
  row.names(aod)[row]
}


.step_p_selected_term <- function(step_label) {
  if (is.na(step_label) || !nzchar(step_label)) {
    return(NA_character_)
  }
  sub("^[+-][[:space:]]*", "", step_label)
}


.step_p_is_cycle <- function(current_label, previous_label) {
  if (is.na(current_label) || is.na(previous_label)) {
    return(FALSE)
  }
  identical(
    .step_p_selected_term(current_label),
    .step_p_selected_term(previous_label)
  )
}



.step_p_backward_candidate <- function(fitted_model,
                                       fscope,
                                       test,
                                       cov_type,
                                       pmethod,
                                       pvalue) {
  if (!length(fscope$drop)) {
    return(NULL)
  }

  aod <- drop1(fitted_model, fscope$drop, test, cov_type, pmethod)
  rn <- row.names(aod)
  row.names(aod) <- c(rn[1L], paste("-", rn[-1L]))
  attr(aod, "heading") <- NULL

  pick_row <- .step_p_selected_row_backward(aod)
  if (.step_p_backward_should_stop(aod, pick_row, pvalue)) {
    return(NULL)
  }

  list(
    aod = aod,
    row = pick_row,
    change = .step_p_selected_step_label(aod, pick_row)
  )
}


.step_p_forward_candidate <- function(fitted_model,
                                      fscope,
                                      test,
                                      cov_type,
                                      pmethod,
                                      pvalue) {
  if (!length(fscope$add)) {
    return(NULL)
  }

  aod <- add1(fitted_model, fscope$add, test, cov_type, pmethod)
  rn <- row.names(aod)
  row.names(aod) <- c(rn[1L], paste("+", rn[-1L]))
  attr(aod, "heading") <- NULL

  pick_row <- .step_p_selected_row_forward(aod)
  if (.step_p_forward_should_stop(aod, pick_row, pvalue)) {
    return(NULL)
  }

  list(
    aod = aod,
    row = pick_row,
    change = .step_p_selected_step_label(aod, pick_row)
  )
}


## backward elimination
step_p_backward <- function(object,
                            scope,
                            test,
                            cov_type,
                            pmethod,
                            pvalue,
                            steps) {
  model_terms <- stats::terms(object)
  object$call$formula <- stats::formula(object)
  scope_factors <- .step_p_scope_factors(object, scope, model_terms)
  factors_drop <- scope_factors$drop
  factors_add <- scope_factors$add
  init <- .step_p_init_models(object, cov_type, steps)
  steps <- init$steps
  models <- init$models
  models_no <- init$models_no
  obs_no <- object$obs_no
  fitted_model <- object
  while (steps > 0L) {
    steps <- steps - 1L
    ffac <- attr(model_terms, "factors")
    fscope <- stats::factor.scope(ffac, list(add = factors_add, drop = factors_drop))
    if (!length(fscope$drop)) {
      break
    }
    aod <- drop1(fitted_model, fscope$drop, test, cov_type, pmethod)
    rn <- row.names(aod)
    row.names(aod) <- c(rn[1L], paste("-", rn[-1L]))
    attr(aod, "heading") <- NULL
    pick_row <- .step_p_selected_row_backward(aod)
    if (.step_p_backward_should_stop(aod, pick_row, pvalue)) {
      break
    }
    change <- row.names(aod)[pick_row]
    fitted_model <- .step_p_update_fit(fitted_model, change, obs_no)
    model_terms <- stats::terms(fitted_model)
    out <- .step_p_append_model(
      models = models,
      models_no = models_no,
      step = change,
      df = aod[pick_row, "Df"],
      chi = aod[pick_row, "Chi"],
      pval = aod[pick_row, "Pr(>Chi)"],
      cic = extract_cic(fitted_model, cov_type)
    )
    models <- out$models
    models_no <- out$models_no
  }
  .step_p_results(models = models[seq_len(models_no)], fit = fitted_model, object = object)
}


## forward selection
step_p_forward <- function(object,
                           scope,
                           test,
                           cov_type,
                           pmethod,
                           pvalue,
                           steps) {
  model_terms <- stats::terms(object)
  object$call$formula <- stats::formula(object)
  scope_factors <- .step_p_scope_factors(object, scope, model_terms)
  factors_drop <- scope_factors$drop
  factors_add <- scope_factors$add
  init <- .step_p_init_models(object, cov_type, steps)
  steps <- init$steps
  models <- init$models
  models_no <- init$models_no
  obs_no <- object$obs_no
  fitted_model <- object

  while (steps > 0L) {
    steps <- steps - 1L
    ffac <- attr(model_terms, "factors")
    fscope <- stats::factor.scope(ffac, list(add = factors_add, drop = factors_drop))
    if (!length(fscope$add)) {
      break
    }
    aod <- add1(fitted_model, fscope$add, test, cov_type, pmethod)
    rn <- row.names(aod)
    row.names(aod) <- c(rn[1L], paste("+", rn[-1L]))
    attr(aod, "heading") <- NULL
    pick_row <- .step_p_selected_row_forward(aod)
    if (.step_p_forward_should_stop(aod, pick_row, pvalue)) {
      break
    }
    change <- row.names(aod)[pick_row]
    fitted_model <- .step_p_update_fit(fitted_model, change, obs_no)
    model_terms <- stats::terms(fitted_model)
    out <- .step_p_append_model(
      models = models,
      models_no = models_no,
      step = change,
      df = aod[pick_row, "Df"],
      chi = aod[pick_row, "Chi"],
      pval = aod[pick_row, "Pr(>Chi)"],
      cic = extract_cic(fitted_model, cov_type)
    )
    models <- out$models
    models_no <- out$models_no
  }

  .step_p_results(models = models[seq_len(models_no)], fit = fitted_model, object = object)
}


## bidirectional elimination
step_p_both <- function(object,
                        scope,
                        test,
                        cov_type,
                        pmethod,
                        p_enter,
                        p_remove,
                        steps) {
  model_terms <- stats::terms(object)
  object$call$formula <- stats::formula(object)

  scope_factors <- .step_p_scope_factors(object, scope, model_terms)
  factors_drop <- scope_factors$drop
  factors_add <- scope_factors$add

  init <- .step_p_init_models(object, cov_type, steps)
  steps <- init$steps
  models <- init$models
  models_no <- init$models_no

  obs_no <- object$obs_no
  fitted_model <- object
  previous_change <- NA_character_

  while (steps > 0L) {
    steps <- steps - 1L

    ffac <- attr(model_terms, "factors")
    fscope <- stats::factor.scope(ffac, list(add = factors_add, drop = factors_drop))

    back <- .step_p_backward_candidate(
      fitted_model = fitted_model,
      fscope = fscope,
      test = test,
      cov_type = cov_type,
      pmethod = pmethod,
      pvalue = p_remove
    )

    candidate <- back

    if (is.null(candidate)) {
      fwd <- .step_p_forward_candidate(
        fitted_model = fitted_model,
        fscope = fscope,
        test = test,
        cov_type = cov_type,
        pmethod = pmethod,
        pvalue = p_enter
      )
      candidate <- fwd
    }

    if (is.null(candidate)) {
      break
    }

    if (.step_p_is_cycle(candidate$change, previous_change)) {
      break
    }

    fitted_model <- .step_p_update_fit(fitted_model, candidate$change, obs_no)
    model_terms <- stats::terms(fitted_model)

    out <- .step_p_append_model(
      models = models,
      models_no = models_no,
      step = candidate$change,
      df = candidate$aod[candidate$row, "Df"],
      chi = candidate$aod[candidate$row, "Chi"],
      pval = candidate$aod[candidate$row, "Pr(>Chi)"],
      cic = extract_cic(fitted_model, cov_type)
    )
    models <- out$models
    models_no <- out$models_no
    previous_change <- candidate$change
  }

  .step_p_results(models = models[seq_len(models_no)], fit = fitted_model, object = object)
}
