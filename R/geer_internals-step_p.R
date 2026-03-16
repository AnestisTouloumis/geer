## backward elimination
step_p_backward <- function(object,
                            scope,
                            test,
                            cov_type,
                            pmethod,
                            pvalue,
                            steps) {
  model_terms <- terms(object)
  object$call$formula <- object$formula <- model_terms
  if (missing(scope)) {
    factors_drop <- numeric()
    factors_add <- attr(model_terms, "factors")
  } else if (is.list(scope)) {
    factors_drop <- scope$lower
    factors_add <- scope$upper
    factors_drop <- if (!is.null(factors_drop)) {
      attr(terms(update.formula(object, factors_drop)), "factors")
    } else {
      numeric()
    }
    factors_add <- if (!is.null(factors_add)) {
      attr(terms(update.formula(object, factors_add)), "factors")
    } else {
      numeric()
    }
  } else {
    factors_add <- if (!is.null(scope)) {
      attr(terms(update.formula(object, scope)), "factors")
    } else {
      numeric()
    }
    factors_drop <- numeric()
  }
  steps <- as.integer(steps)
  if (length(steps) != 1L || steps < 0L) {
    stop("'steps' must be a non-negative integer", call. = FALSE)
  }
  models_no <- 1L
  models <- vector("list", steps + 1L)
  models[[models_no]] <- list(
    Step = "",
    Df = NA_real_,
    Chi = NA_real_,
    `Pr(>Chi)` = NA_real_,
    CIC = extract_cic(object, cov_type)
  )
  obs_no <- object$obs_no
  fitted_model <- object
  while (steps > 0L) {
    steps <- steps - 1L
    ffac <- attr(model_terms, "factors")
    fscope <- factor.scope(ffac, list(add = factors_add, drop = factors_drop))
    if (!length(fscope$drop)) break
    aod <- drop1(fitted_model, fscope$drop, test, cov_type, pmethod)
    rn <- row.names(aod)
    row.names(aod) <- c(rn[1L], paste("-" , rn[-1L]))
    attr(aod, "heading") <- NULL
    aod2 <- aod[-1L, , drop = FALSE]
    if (!nrow(aod2)) break
    zdf <- !is.na(aod2$Df) & aod2$Df == 0
    if (any(zdf)) {
      change <- row.names(aod2)[which(zdf)][length(which(zdf))]
      pick_row <- match(change, row.names(aod))
    } else {
      pvals <- aod2[, "Pr(>Chi)"]
      ok <- is.finite(pvals)
      if (!any(ok)) break
      idx <- which.max(pvals[ok])
      change <- row.names(aod2)[which(ok)[idx]]
      if (pvals[which(ok)[idx]] < pvalue) break
      pick_row <- match(change, row.names(aod))
    }
    fitted_model <- update(
      fitted_model,
      formula = paste(". ~ .", change),
      evaluate = FALSE
    )
    fitted_model <- eval.parent(fitted_model)
    obs_no_new <- fitted_model$obs_no
    if (all(is.finite(c(obs_no, obs_no_new))) && obs_no_new != obs_no) {
      stop("number of rows in use has changed: remove missing values?", call. = FALSE)
    }
    model_terms <- terms(fitted_model)
    models_no <- models_no + 1L
    models[[models_no]] <- list(
      Step = change,
      Df = aod[pick_row, "Df"],
      Chi = aod[pick_row, "Chi"],
      `Pr(>Chi)` = aod[pick_row, "Pr(>Chi)"],
      CIC = extract_cic(fitted_model, cov_type)
    )
  }
  step_results(models = models[seq_len(models_no)], fit = fitted_model, object = object)
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
  object$call$formula <- object$formula <- model_terms
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
  steps <- as.integer(steps)
  if (length(steps) != 1L || steps < 0L) {
    stop("'steps' must be a non-negative integer", call. = FALSE)
  }
  models_no <- 1L
  models <- vector("list", steps + 1L)
  models[[models_no]] <- list(
    Step = "",
    Df = NA_real_,
    Chi = NA_real_,
    `Pr(>Chi)` = NA_real_,
    CIC = extract_cic(object, cov_type)
  )
  obs_no <- object$obs_no
  fitted_model <- object
  while (steps > 0L) {
    steps <- steps - 1L
    ffac <- attr(model_terms, "factors")
    fscope <- factor.scope(ffac, list(add = factors_add, drop = factors_drop))
    if (!length(fscope$add)) break
    aod <- add1(fitted_model, fscope$add, test, cov_type, pmethod)
    rn <- row.names(aod)
    row.names(aod) <- c(rn[1L], paste("+", rn[-1L]))
    attr(aod, "heading") <- NULL
    aod2 <- aod[-1L, , drop = FALSE]
    if (!nrow(aod2)) break
    keep <- is.na(aod2$Df) | aod2$Df != 0
    aod2 <- aod2[keep, , drop = FALSE]
    if (!nrow(aod2)) break
    pvals <- aod2[, "Pr(>Chi)"]
    ok <- is.finite(pvals)
    if (!any(ok)) break
    idx <- which.min(pvals[ok])
    change <- row.names(aod2)[which(ok)[idx]]
    if (pvals[which(ok)[idx]] > pvalue) break
    pick_row <- match(change, row.names(aod))
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
    model_terms <- stats::terms(fitted_model)
    models_no <- models_no + 1L
    models[[models_no]] <- list(
      Step = change,
      Df = aod[pick_row, "Df"],
      Chi = aod[pick_row, "Chi"],
      `Pr(>Chi)` = aod[pick_row, "Pr(>Chi)"],
      CIC = extract_cic(fitted_model, cov_type)
    )
  }
  step_results(models = models[seq_len(models_no)], fit = fitted_model, object = object)
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
  steps <- as.integer(steps)
  if (length(steps) != 1L || steps < 0L) {
    stop("'steps' must be a non-negative integer", call. = FALSE)
  }
  models_no <- 1L
  models <- vector("list", steps + 1L)
  models[[models_no]] <- list(
    Step = "",
    Df = NA_real_,
    Chi = NA_real_,
    `Pr(>Chi)` = NA_real_,
    CIC = extract_cic(object, cov_type)
  )
  fitted_model <- object
  same_variable <- 0L
  last_change <- ""
  while (steps > 0L && same_variable < 2L) {
    ## backward (one step)
    fitted_model <- step_p_backward(
      fitted_model,
      scope = scope,
      test = test,
      cov_type = cov_type,
      pmethod = pmethod,
      pvalue = p_remove,
      steps = 1
    )
    steps <- steps - 1L
    moved <- !all(is.na(fitted_model$anova$`Pr(>Chi)`))
    if (!moved) {
      same_variable <- same_variable + 1L
    } else {
      same_variable <- 0L
      models_no <- models_no + 1L
      models[[models_no]] <- list(
        Step = fitted_model$anova$Step[2],
        Df = fitted_model$anova$Df[2],
        Chi = fitted_model$anova$Chi[2],
        `Pr(>Chi)` = fitted_model$anova$`Pr(>Chi)`[2],
        CIC = fitted_model$anova$CIC[2]
      )
      cur <- sub("^-\\s*", "", models[[models_no]]$Step)
      if (identical(last_change, cur)) break
      last_change <- cur
    }
    if (steps <= 0L) break
    fitted_model <- step_p_forward(
      fitted_model,
      scope = scope,
      test = test,
      cov_type = cov_type,
      pmethod = pmethod,
      pvalue = p_enter,
      steps = 1
    )
    steps <- steps - 1L
    moved <- !all(is.na(fitted_model$anova$`Pr(>Chi)`))
    if (!moved) {
      same_variable <- same_variable + 1L
    } else {
      same_variable <- 0L
      models_no <- models_no + 1L
      models[[models_no]] <- list(
        Step = fitted_model$anova$Step[2],
        Df = fitted_model$anova$Df[2],
        Chi = fitted_model$anova$Chi[2],
        `Pr(>Chi)` = fitted_model$anova$`Pr(>Chi)`[2],
        CIC = fitted_model$anova$CIC[2]
      )
      cur <- sub("^\\+\\s*", "", models[[models_no]]$Step)
      if (identical(last_change, cur)) break
      last_change <- cur
    }
  }
  step_results(models = models[seq_len(models_no)], fit = fitted_model, object = object)
}


## printing results
step_results <- function(models, fit, object) {
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
