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
  } else {
    if (is.list(scope)) {
      factors_drop <- scope$lower
      factors_add <- scope$upper
      factors_drop <- if (!is.null(factors_drop)) {
        attr(terms(update.formula(object, factors_drop)), "factors")
      } else {
        numeric()
      }
      if (!is.null(factors_add))
        factors_add <- attr(terms(update.formula(object, factors_add)), "factors")
    } else {
      factors_add <- scope
      if (!is.null(factors_add))
        factors_add <- attr(terms(update.formula(object, scope)), "factors")
      factors_drop <- numeric()
    }
  }
  models_no <- 1
  models <- vector("list", steps)
  models[[models_no]] <-
    list(Step = "", Df = NA, Chi = NA, `Pr(>Chi)` = NA,
         CIC = extract_cic(object, cov_type))
  obs_no <- object$obs_no
  fitted_model <- object
  while (steps > 0) {
    steps <- steps - 1
    ffac <- attr(model_terms, "factors")
    scope <- factor.scope(ffac, list(add = factors_add, drop = factors_drop))
    aod <- change <- NULL
    if (length(scope$drop)) {
      aod <- drop1(fitted_model, scope$drop, test, cov_type, pmethod)
      rn <- row.names(aod)
      row.names(aod) <- c(rn[1L], paste("-", rn[-1L]))
      if (any(aod$Df == 0, na.rm = TRUE)) {
        zdf <- aod$Df == 0 & !is.na(aod$Df)
        change <- rev(rownames(aod)[zdf])[1L]
      }
    }
    if (is.null(change)) {
      attr(aod, "heading") <- NULL
      nzdf <- if (!is.null(aod$Df))
        aod$Df != 0 | !is.na(aod$Df)
      aod <- aod[nzdf, ]
      if (is.null(aod) || ncol(aod) == 0)
        break
      pvalue_max <- which.max(aod[, "Pr(>Chi)"])
      if (pvalue_max == 1)
        break
      criterion <- aod[pvalue_max, "Pr(>Chi)"] < pvalue
      if (criterion)
        break
      change <- rownames(aod)[pvalue_max]
    }
    fitted_model <- update(fitted_model,
                           formula = paste(". ~ .", change),
                           evaluate = FALSE)
    fitted_model <- eval.parent(fitted_model)
    obs_no_new <- fitted_model$obs_no
    if (all(is.finite(c(obs_no, obs_no_new))) && obs_no_new != obs_no)
      stop("number of rows in use has changed: remove missing values?")
    model_terms <- terms(fitted_model)
    models_no <- models_no + 1
    models[[models_no]] <-
      list(Step = change, Df = aod[pvalue_max, "Df"], Chi = aod[pvalue_max, "Chi"],
           `Pr(>Chi)` = aod[pvalue_max, "Pr(>Chi)"],
           CIC = extract_cic(fitted_model, cov_type))
  }
  step_results(models = models[seq(models_no)], fitted_model, object)
}


## forward selection
step_p_forward <- function(object,
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
  } else {
    if (is.list(scope)) {
      factors_drop <- scope$lower
      factors_add <- scope$upper
      factors_drop <- if (!is.null(factors_drop)) {
        attr(terms(update.formula(object, factors_drop)), "factors")
      } else {
        numeric()
      }
      if (!is.null(factors_add))
        factors_add <- attr(terms(update.formula(object, factors_add)), "factors")
    } else {
      factors_add <- scope
      if (!is.null(factors_add))
        factors_add <- attr(terms(update.formula(object, scope)), "factors")
      factors_drop <- numeric()
    }
  }
  models_no <- 1
  models <- vector("list", steps)
  models[[models_no]] <-
    list(Step = "", Df = NA, Chi = NA, `Pr(>Chi)` = NA,
         CIC = extract_cic(object, cov_type))
  obs_no <- object$obs_no
  fitted_model <- object
  while (steps > 0) {
    steps <- steps - 1
    ffac <- attr(model_terms, "factors")
    scope <- factor.scope(ffac, list(add = factors_add, drop = factors_drop))
    change <- NULL
    if (length(scope$add)) {
      aod  <- add1(fitted_model, scope$add, test, cov_type, pmethod)
      rn <- row.names(aod)
      row.names(aod) <- c(rn[1L], paste("+", rn[-1L]))
    } else break
    attr(aod, "heading") <- NULL
    nzdf <- if (!is.null(aod$Df))
      aod$Df != 0 | is.na(aod$Df)
    aod <- aod[nzdf, ]
    if (is.null(aod) || ncol(aod) == 0)
      break
    pvalue_min <- which.min(aod[, "Pr(>Chi)"])
    if (pvalue_min == 1)
      break
    criterion <- aod[pvalue_min, "Pr(>Chi)"] > pvalue
    if (criterion)
      break
    change <- rownames(aod)[pvalue_min]
    fitted_model <- update(fitted_model,
                           formula = paste(" . ~ .", change),
                           evaluate = FALSE)
    fitted_model <- eval.parent(fitted_model)
    obs_no_new <- fitted_model$obs_no
    if (all(is.finite(c(obs_no, obs_no_new))) && obs_no_new != obs_no)
      stop("number of rows in use has changed: remove missing values?")
    model_terms <- terms(fitted_model)
    models_no <- models_no + 1
    models[[models_no]] <-
      list(Step = change, Df = aod[pvalue_min, "Df"],
           Chi = aod[pvalue_min, "Chi"],
           `Pr(>Chi)` = aod[pvalue_min, "Pr(>Chi)"],
           CIC = extract_cic(fitted_model, cov_type))
  }
  step_results(models = models[seq(models_no)], fitted_model, object)
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
  models <- vector("list", steps)
  models_no <- 1
  models <- vector("list", steps)
  models[[models_no]] <-
    list(Step = "", Df = NA, Chi = NA, `Pr(>Chi)` = NA,
         CIC = extract_cic(object, cov_type))
  fitted_model <- object
  same_variable <- 0
  change_variable <- ""
  while (steps > 0 & same_variable < 2) {
    steps <- steps - 1
    fitted_model <- step_p_backward(fitted_model,
                                    scope = scope,
                                    test = test,
                                    cov_type = cov_type,
                                    pmethod = pmethod,
                                    pvalue = p_remove,
                                    steps = 1)
    if (all(is.na(fitted_model$anova$`Pr(>Chi)`))) {
      same_variable <- same_variable + 1
    } else {
      same_variable <- 0
      models_no <- models_no + 1
      models[[models_no]] <-
        list(Step = fitted_model$anova$Step[2],
             Df = fitted_model$anova$Df[2],
             Chi = fitted_model$anova$Chi[2],
             `Pr(>Chi)` = fitted_model$anova$`Pr(>Chi)`[2],
             CIC = fitted_model$anova$CIC[2])
      if (change_variable == sub("- ", "", models[[models_no]]$Step)) {
        break
      } else {
        change_variable <- sub("- ", "", models[[models_no]]$Step)
      }
    }
    steps <- steps - 1
    fitted_model <- step_p_forward(fitted_model,
                                   scope = scope,
                                   test = test,
                                   cov_type = cov_type,
                                   pmethod = pmethod,
                                   pvalue = p_enter,
                                   steps = 1)
    if (all(is.na(fitted_model$anova$`Pr(>Chi)`))) {
      same_variable <- same_variable + 1
    } else {
      same_variable <- 0
      models_no <- models_no + 1
      models[[models_no]] <-
        list(Step = fitted_model$anova$Step[2],
             Df = fitted_model$anova$Df[2],
             Chi = fitted_model$anova$Chi[2],
             `Pr(>Chi)` = fitted_model$anova$`Pr(>Chi)`[2],
             CIC = fitted_model$anova$CIC[2])
      if (change_variable == sub("+ ", "", models[[models_no]]$Step)) {
        break
      } else {
        change_variable <- sub("+ ", "", models[[models_no]]$Step)
      }
    }
  }
  step_results(models = models[seq(models_no)], fitted_model, object)
}


## printing results
step_results <- function(models, fit, object) {
  Step <- sapply(models, `[[`, "Step")
  Df <- sapply(models, `[[`, "Df")
  Chi <- sapply(models, `[[`, "Chi")
  pvalue <- sapply(models, `[[`, "Pr(>Chi)")
  CIC <- sapply(models, `[[`, "CIC")
  heading <- c("Stepwise Model Path \nAnalysis of Deviance Table",
               "\nInitial Model:", deparse(formula(object)),
               "\nFinal Model:", deparse(formula(fit)), "\n")
  aod <- data.frame(Step = Step, Df = Df, Chi = Chi, `Pr(>Chi)` = pvalue,
                    CIC = CIC, check.names = FALSE)
  attr(aod, "heading") <- heading
  structure(aod, title = heading, class = c("anova", "data.frame"))
  fit$anova <- aod
  fit
}


