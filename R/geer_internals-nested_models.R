## checking if two models are nested
check_nested_models <- function(object0, object1) {
  if ( !("geer" %in% class(object0)) | !("geer" %in% class(object1)) ) {
    stop("Objects must be of 'geer' class", call. = FALSE)
  }
  if (!all(object0$y == object1$y)) {
    stop("Response variable differs in the models", call. = FALSE)
  }
  n0 <- length(object0$coefficients)
  n1 <- length(object1$coefficients)
  if (n0 == n1) {
    stop("Models must be nested", call. = FALSE)
  }
  if (n0 < n1) {
    obj0 <- object0
    obj1 <- object1
  }
  else {
    obj0 <- object1
    obj1 <- object0
  }
  names0 <- names(obj0$coefficients)
  names1 <- names(obj1$coefficients)
  names_test <- setdiff(names1, names0)
  if (length(names_test) == 0) {
    stop("Models must be nested", call. = FALSE)
  }
  if (length(setdiff(names0, names1)) != 0) {
    stop("Models must be nested", call. = FALSE)
  }
  index <- rep(0, length(names_test))
  for (i in seq_len(length(names_test))) {
    index[i] <- which(names_test[i] == names1)
  }
  list(obj0 = obj0, obj1 = obj1, index = index)
}


## chisq-approximations for calculating the p-value for sum of independent
## chi sq random variables
lcsumchisq <- function(x, test_stat, pmethod) {
  x_bar <- mean(x)
  df <- length(x)
  if (pmethod == "rao-scott") {
    test_stat <- test_stat/x_bar
    pvalue <- 1 - pchisq(test_stat, df = df)
  } else if (pmethod == "satterthwaite") {
    alpha <- sum((x - x_bar)^2)/(df * x_bar^2)
    df <- df/(1 + alpha^2)
    test_stat <- test_stat/((1 + alpha^2) * x_bar)
    pvalue <- 1 - pchisq(test_stat, df = df)
  }
  list(test_stat = test_stat, df = df, pvalue = pvalue)
}


## wald test
wald_test <- function(object0, object1, cov_type = "robust"){
  nested_models <- check_nested_models(object0, object1)
  obj0 <- nested_models$obj0
  obj1 <- nested_models$obj1
  index <- nested_models$index
  df <- length(index)
  coeffs_test <- obj1$coefficients[index]
  cov_test <- vcov(nested_models$obj1, cov_type = cov_type)[index, index]
  wald_stat <- t(coeffs_test) %*% solve(cov_test, coeffs_test)
  pvalue <- 1 - pchisq(wald_stat, df)
  ans <- list(Df = df, X2 = wald_stat, `P(>X2)` = pvalue)
  ans
}


## working wald test
working_wald_test <- function(object0,
                              object1,
                              cov_type = "robust",
                              pmethod = "rao-scott"){
  nested_models <- check_nested_models(object0, object1)
  obj1 <- nested_models$obj1
  obj0 <- nested_models$obj0
  index <- nested_models$index
  coeffs_test <- obj1$coefficients[index]
  cov_test <- vcov(obj1, cov_type = "naive")[index, index]
  working_wald_stat <- t(coeffs_test) %*% solve(cov_test, coeffs_test)
  rob_test <- vcov(obj1, cov_type = cov_type)[index, index]
  coeffs <- eigen(solve(cov_test, rob_test), only.values = TRUE)$values
  ans_stats <- lcsumchisq(coeffs, working_wald_stat, pmethod)
  ans <- list(Df = ans_stats$df, X2 = ans_stats$test_stat, `P(>X2)` = ans_stats$pvalue)
  ans
}


## working likelihood ratio test
working_lr_test <- function(object0,
                            object1,
                            cov_type = "robust",
                            pmethod = "rao-scott"){
  nested_models <- check_nested_models(object0, object1)
  obj1 <- nested_models$obj1
  obj0 <- nested_models$obj0
  index <- nested_models$index
  ll_obj0 <- compute_quasi_loglikelihood(obj0)
  ll_obj1 <- compute_quasi_loglikelihood(obj1)
  if (!all(c(obj0$association_structure, obj1$association_structure) == "independence"))
    stop("the working-LRT is available only for independence working models")
  working_ll_test <- -2 * (ll_obj0 - ll_obj1)
  cov_test <- vcov(obj1, cov_type = "naive")[index, index]
  rob_test <- vcov(obj1, cov_type = cov_type)[index, index]
  coeffs <- eigen(solve(cov_test, rob_test), only.values = TRUE)$values
  ans_stats <- lcsumchisq(coeffs, working_ll_test, pmethod)
  ans <- list(Df = ans_stats$df, X2 = ans_stats$test_stat, `P(>X2)` = ans_stats$pvalue)
  ans
}

## score test
score_test <- function(object0, object1, cov_type = "robust"){
  nested_models <- check_nested_models(object0, object1)
  obj1 <- nested_models$obj1
  obj0 <- nested_models$obj0
  index <- nested_models$index
  df <- length(index)
  coeffs_test <- obj1$coefficients
  coeffs_test[index] <- 0
  coeffs_test[-index] <- obj0$coefficients
  if (obj1$call[1] == "geewa()") {
    uvector <- estimating_equations_gee_cc(obj1$y,
                                           obj1$model_matrix,
                                           obj1$id,
                                           obj1$repeated,
                                           obj1$weights,
                                           obj1$family$link,
                                           obj1$family$family,
                                           coeffs_test,
                                           obj0$fitted.values,
                                           obj0$linear.predictors,
                                           obj1$association_structure,
                                           obj1$alpha,
                                           obj1$phi)

    covariance <- get_covariance_matrices_cc(obj1$y,
                                             obj1$model_matrix,
                                             obj1$id,
                                             obj1$repeated,
                                             obj1$weights,
                                             obj1$family$link,
                                             obj1$family$family,
                                             obj0$fitted.values,
                                             obj0$linear.predictors,
                                             obj1$association_structure,
                                             obj1$alpha,
                                             obj1$phi)
  } else {
    if (length(obj1$alpha) == 1) {
      association_alpha <- rep(obj1$alpha,
                               choose(max(obj1$repeated), 2))
    }
    uvector <- estimating_equations_gee_or(obj1$y,
                                           obj1$model_matrix,
                                           obj1$id,
                                           obj1$repeated,
                                           obj1$weights,
                                           obj1$family$link,
                                           coeffs_test,
                                           obj0$fitted.values,
                                           obj0$linear.predictors,
                                           association_alpha)
    covariance <- get_covariance_matrices_or(obj1$y,
                                             obj1$model_matrix,
                                             obj1$id,
                                             obj1$repeated,
                                             obj0$weights,
                                             obj1$family$link,
                                             obj0$fitted.values,
                                             obj0$linear.predictors,
                                             association_alpha)
  }
  if (cov_type == "robust") {
    cov_test <- covariance$robust_covariance[index, index]
  } else if (cov_type == "naive") {
    cov_test <- covariance$naive_covariance[index, index]
  } else if (cov_type == "bias-corrected") {
    cov_test <- covariance$bc_covariance[index, index]
  } else {
    sample_size <- obj1$clusters_no
    parameters_no <- length(coef(obj1))
    cov_test <-
      (sample_size / (sample_size - parameters_no)) *
      covariance$robust_covariance[index, index]
  }
  naive_cov <-  covariance$naive_covariance
  score_stat <-
    (t(uvector) %*% naive_cov[, index]) %*%
    (solve(cov_test) %*% (naive_cov[index, ] %*% uvector))
  pvalue <- 1 - pchisq(score_stat, df)
  ans <- list(Df = df, X2 = score_stat, `P(>X2)` = pvalue)
  ans
}

## working score test
working_score_test <- function(object0, object1, cov_type = "robust", pmethod = "rao-scott"){
  nested_models <- check_nested_models(object0, object1)
  obj1 <- nested_models$obj1
  obj0 <- nested_models$obj0
  index <- nested_models$index
  coeffs_test <- obj1$coefficients
  coeffs_test[index] <- 0
  coeffs_test[-index] <- obj0$coefficients
  if (obj1$call[1] == "geewa()") {
    uvector <- estimating_equations_gee_cc(obj1$y,
                                           obj1$model_matrix,
                                           obj1$id,
                                           obj1$repeated,
                                           obj1$weights,
                                           obj1$family$link,
                                           obj1$family$family,
                                           coeffs_test,
                                           obj0$fitted.values,
                                           obj0$linear.predictors,
                                           obj1$association_structure,
                                           obj1$alpha,
                                           obj1$phi)
    covariance <- get_covariance_matrices_cc(obj1$y,
                                             obj1$model_matrix,
                                             obj1$id,
                                             obj1$repeated,
                                             obj1$weights,
                                             obj1$family$link,
                                             obj1$family$family,
                                             obj0$fitted.values,
                                             obj0$linear.predictors,
                                             obj1$association_structure,
                                             obj1$alpha,
                                             obj1$phi)
  } else {
    if (length(obj1$alpha) == 1) {
      association_alpha <- rep(obj1$alpha,
                               choose(max(obj1$repeated), 2))
    }
    uvector <- estimating_equations_gee_or(obj1$y,
                                           obj1$model_matrix,
                                           obj1$id,
                                           obj1$repeated,
                                           obj1$weights,
                                           obj1$family$link,
                                           coeffs_test,
                                           obj0$fitted.values,
                                           obj0$linear.predictors,
                                           association_alpha)
    covariance <- get_covariance_matrices_or(obj1$y,
                                             obj1$model_matrix,
                                             obj1$id,
                                             obj1$repeated,
                                             obj0$weights,
                                             obj1$family$link,
                                             obj0$fitted.values,
                                             obj0$linear.predictors,
                                             association_alpha)
  }
  cov_test <- covariance$naive_covariance
  if (cov_type == "robust") {
    rob_test <- covariance$robust_covariance
  } else if (cov_type == "naive") {
    rob_test <- covariance$naive_covariance
  } else if (cov_type == "bias-corrected") {
    rob_test <- covariance$bc_covariance
  } else {
    sample_size <- obj1$clusters_no
    parameters_no <- length(coef(obj1))
    rob_test <-
      (sample_size / (sample_size - parameters_no)) * covariance$robust_covariance
  }
  score_stat <-
    (t(uvector) %*% cov_test[, index]) %*%
    (solve(cov_test[index, index]) %*% (cov_test[index, ] %*% uvector))
  coeffs <-
    eigen(solve(cov_test[index, index], rob_test[index, index]), only.values = TRUE)$values
  ans_stats <- lcsumchisq(coeffs, score_stat, pmethod)
  ans <- list(Df = ans_stats$df, X2 = ans_stats$test_stat, `P(>X2)` = ans_stats$pvalue)
  ans
}

anova.geerlist <-
  function(object, ...,
           test = "wald",
           cov_type = "robust",
           pmethod = "rao-scott")
  {
    response_vectors <-
      as.character(lapply(object, function(x) deparse(formula(x)[[2L]])))
    same_response <- response_vectors == response_vectors[1L]
    if (!all(same_response)) {
      object <- object[same_response]
      warning("Models with response ", deparse(response_vectors[!same_response]),
              " removed because response differs from Model 1")
    }
    association_vector <-
      as.character(lapply(object, function(x) x$association_structure))
    same_association <- association_vector == association_vector[1L]
    if (!all(same_association)) {
      object <- object[same_association]
      warning("Models with association structure ",
              deparse(association_vector[!same_association]),
              " removed because association structure differs from ", "model 1")
    }
    sample_size_vector <- sapply(object, function(x) length(x$residuals))
    if (any(sample_size_vector != sample_size_vector[1L]))
      stop("models were not all fitted to the same size of dataset")
    models_no <- length(object)
    if (models_no == 1)
      return(anova.geer(object[[1L]], cov_type = cov_type, pmethod = pmethod))
    variables <- lapply(object, function(x)
      paste(deparse(formula(x)), collapse = "\n") )
    table <- matrix(0, models_no - 1, 3)
    test_type <- switch(test,
                        wald = "Wald",
                        score = "Score",
                        `working-wald` = "Working Wald",
                        `working-score` = "Working Score",
                        `working-lrt` = "Working LRT")
    for (i in 1:nrow(table)) {
      value <-
        switch(test,
               wald =
                 wald_test(object[[i]], object[[i + 1]], cov_type),
               score =
                 score_test(object[[i]], object[[i + 1]], cov_type),
               `working-wald` =
                 working_wald_test(object[[i]], object[[i + 1]], cov_type, pmethod),
               `working-score` =
                 working_score_test(object[[i]], object[[i + 1]], cov_type, pmethod),
               `working-lrt` =
                 working_lr_test(object[[i]], object[[i + 1]], cov_type, pmethod)
        )
      table[i, ] <- c(value$Df, value$X2, value$`P(>X2)`)
    }
    title  <- paste("Analysis of", test_type, "Statistic Table\n", sep = " ")
    topnote <-
      paste("Model ",
            format(1L:models_no),
            ": ",
            variables,
            sep = "",
            collapse = "\n")
    table <- data.frame(table)
    dimnames(table) <- list(paste(1L:(models_no - 1), "vs", 2L:models_no),
                            c("Df", "X2", "P(>X2)"))
    structure(table,
              heading = c(title, topnote),
              class = c("anova", "data.frame"))
  }
