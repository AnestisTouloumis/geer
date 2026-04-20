.step_p_change_label <- function(term, action = c("add", "drop")) {
  action <- match.arg(action)
  prefix <- if (identical(action, "add")) "+ " else "- "
  paste0(prefix, term)
}


.step_p_validate_scope <- function(scope) {
  if (is.null(scope) || inherits(scope, "formula") || is.character(scope)) {
    return(scope)
  }
  if (!is.list(scope)) {
    stop("'scope' must be NULL, a formula, a character vector, or a list",
         call. = FALSE)
  }
  scope_names <- names(scope)
  if (is.null(scope_names) || !any(scope_names %in% c("lower", "upper"))) {
    stop("list 'scope' must contain a 'lower' component, an 'upper' component, or both",
         call. = FALSE)
  }
  bad_names <- setdiff(scope_names, c("lower", "upper", ""))
  if (length(bad_names)) {
    stop("list 'scope' can only contain 'lower' and 'upper' components",
         call. = FALSE)
  }
  for (nm in c("lower", "upper")) {
    x <- scope[[nm]]
    if (!is.null(x) && !inherits(x, "formula") && !is.character(x)) {
      stop(sprintf("'scope$%s' must be NULL, a formula, or a character vector", nm),
           call. = FALSE)
    }
  }
  scope
}


.step_p_scope_terms <- function(base_formula, scope_term) {
  if (is.null(scope_term)) {
    return(numeric())
  }
  attr(stats::terms(stats::update.formula(base_formula, scope_term)), "factors")
}


.step_p_scope_factors <- function(object, scope, model_terms) {
  base_formula <- stats::formula(object)
  scope <- .step_p_validate_scope(scope)
  if (is.null(scope)) {
    return(list(
      add = attr(model_terms, "factors"),
      drop = numeric()
    ))
  }
  if (inherits(scope, "formula") || is.character(scope)) {
    return(list(
      add = .step_p_scope_terms(base_formula, scope),
      drop = numeric()
    ))
  }
  list(
    add = .step_p_scope_terms(base_formula, scope$upper),
    drop = .step_p_scope_terms(base_formula, scope$lower)
  )
}


.step_p_init_models <- function(object, cov_type, steps) {
  steps <- as.integer(steps)
  models <- vector("list", steps + 1L)
  models[[1L]] <- list(
    Step = 0L,
    Change = NA_character_,
    Df = NA_real_,
    Chi = NA_real_,
    `Pr(>Chi)` = NA_real_,
    CIC = compute_gee_cic(object, cov_type)
  )
  list(
    steps = steps,
    models = models,
    models_no = 1L
  )
}


.step_p_update_fit <- function(fitted_model, change, obs_no) {
  updated_model <- update(
    fitted_model,
    formula = paste(". ~ .", change),
    data = fitted_model$data
  )
  updated_model <- restore_original_data_call(updated_model, fitted_model)
  obs_no_new <- updated_model$obs_no
  if (all(is.finite(c(obs_no, obs_no_new))) && obs_no_new != obs_no) {
    stop("number of rows in use has changed: remove missing values?", call. = FALSE)
  }
  updated_model
}


.step_p_append_model <- function(models, models_no, change, df, chi, pval, cic) {
  models_no <- models_no + 1L
  models[[models_no]] <- list(
    Step = models_no - 1L,
    Change = change,
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


.step_p_results <- function(models, fit, object) {
  step_no <- vapply(models, `[[`, integer(1), "Step")
  change <- vapply(models, function(x) {
    out <- x[["Change"]]
    if (is.null(out) || is.na(out)) "" else out
  }, character(1))
  row_labels <- ifelse(
    step_no == 0L,
    "Initial model",
    paste0("Step ", step_no, ": ", change)
  )
  aod <- data.frame(
    Step = ifelse(step_no == 0L, "", as.character(step_no)),
    Df = vapply(models, `[[`, NA_real_, "Df"),
    Chi = vapply(models, `[[`, NA_real_, "Chi"),
    `Pr(>Chi)` = vapply(models, `[[`, NA_real_, "Pr(>Chi)"),
    CIC = vapply(models, `[[`, NA_real_, "CIC"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  row.names(aod) <- row_labels
  heading <- c(
    "Stepwise Model Path \nAnalysis of Deviance Table",
    "\nInitial Model:", paste(deparse(stats::formula(object)), collapse = " "),
    "\nFinal Model:", paste(deparse(stats::formula(fit)), collapse = " "), "\n"
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
  aod <- drop1(
    object = fitted_model,
    scope = fscope$drop,
    test = test,
    cov_type = cov_type,
    pmethod = pmethod
  )
  attr(aod, "heading") <- NULL

  pick_row <- .step_p_selected_row_backward(aod)
  if (.step_p_backward_should_stop(aod, pick_row, pvalue)) {
    return(NULL)
  }
  selected_term <- fscope$drop[[pick_row - 1L]]
  list(
    aod = aod,
    row = pick_row,
    change = .step_p_change_label(selected_term, action = "drop")
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
  aod <- add1(
    object = fitted_model,
    scope = fscope$add,
    test = test,
    cov_type = cov_type,
    pmethod = pmethod
  )
  attr(aod, "heading") <- NULL
  pick_row <- .step_p_selected_row_forward(aod)
  if (.step_p_forward_should_stop(aod, pick_row, pvalue)) {
    return(NULL)
  }
  selected_term <- fscope$add[[pick_row - 1L]]
  list(
    aod = aod,
    row = pick_row,
    change = .step_p_change_label(selected_term, action = "add")
  )
}


.step_p_run_backward <- function(object,
                                 scope,
                                 test,
                                 cov_type,
                                 pmethod,
                                 pvalue,
                                 steps) {
  model_terms <- stats::terms(object)
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
    candidate <- .step_p_backward_candidate(
      fitted_model = fitted_model,
      fscope = fscope,
      test = test,
      cov_type = cov_type,
      pmethod = pmethod,
      pvalue = pvalue
    )
    if (is.null(candidate)) {
      break
    }
    fitted_model <- .step_p_update_fit(fitted_model, candidate$change, obs_no)
    model_terms <- stats::terms(fitted_model)
    out <- .step_p_append_model(
      models = models,
      models_no = models_no,
      change = candidate$change,
      df = candidate$aod[candidate$row, "Df"],
      chi = candidate$aod[candidate$row, "Chi"],
      pval = candidate$aod[candidate$row, "Pr(>Chi)"],
      cic = compute_gee_cic(fitted_model, cov_type)
    )
    models <- out$models
    models_no <- out$models_no
  }
  .step_p_results(models = models[seq_len(models_no)], fit = fitted_model, object = object)
}


.step_p_run_forward <- function(object,
                                scope,
                                test,
                                cov_type,
                                pmethod,
                                pvalue,
                                steps) {
  model_terms <- stats::terms(object)
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
    candidate <- .step_p_forward_candidate(
      fitted_model = fitted_model,
      fscope = fscope,
      test = test,
      cov_type = cov_type,
      pmethod = pmethod,
      pvalue = pvalue
    )
    if (is.null(candidate)) {
      break
    }
    fitted_model <- .step_p_update_fit(fitted_model, candidate$change, obs_no)
    model_terms <- stats::terms(fitted_model)
    out <- .step_p_append_model(
      models = models,
      models_no = models_no,
      change = candidate$change,
      df = candidate$aod[candidate$row, "Df"],
      chi = candidate$aod[candidate$row, "Chi"],
      pval = candidate$aod[candidate$row, "Pr(>Chi)"],
      cic = compute_gee_cic(fitted_model, cov_type)
    )
    models <- out$models
    models_no <- out$models_no
  }
  .step_p_results(models = models[seq_len(models_no)], fit = fitted_model, object = object)
}


.step_p_run_both <- function(object,
                             scope,
                             test,
                             cov_type,
                             pmethod,
                             p_enter,
                             p_remove,
                             steps) {
  model_terms <- stats::terms(object)
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
    candidate <- .step_p_backward_candidate(
      fitted_model = fitted_model,
      fscope = fscope,
      test = test,
      cov_type = cov_type,
      pmethod = pmethod,
      pvalue = p_remove
    )
    if (is.null(candidate)) {
      candidate <- .step_p_forward_candidate(
        fitted_model = fitted_model,
        fscope = fscope,
        test = test,
        cov_type = cov_type,
        pmethod = pmethod,
        pvalue = p_enter
      )
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
      change = candidate$change,
      df = candidate$aod[candidate$row, "Df"],
      chi = candidate$aod[candidate$row, "Chi"],
      pval = candidate$aod[candidate$row, "Pr(>Chi)"],
      cic = compute_gee_cic(fitted_model, cov_type)
    )
    models <- out$models
    models_no <- out$models_no
    previous_change <- candidate$change
  }
  .step_p_results(models = models[seq_len(models_no)], fit = fitted_model, object = object)
}
