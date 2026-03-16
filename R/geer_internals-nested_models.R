## checking if two models are nested
check_nested_models <- function(object0, object1) {
  if (!inherits(object0, "geer") || !inherits(object1, "geer")) {
    stop("Objects must be of 'geer' class", call. = FALSE)
  }
  if (!is.null(object0$obs_no) && !is.null(object1$obs_no)) {
    if (!identical(object0$obs_no, object1$obs_no)) {
      stop("Models were not fit on the same observations", call. = FALSE)
    }
  } else {
    if (!identical(object0$y, object1$y)) {
      stop("Response variable differs in the models", call. = FALSE)
    }
    if (!identical(object0$id, object1$id)) {
      stop("Cluster identifiers differs in the models", call. = FALSE)
    }
    if (!identical(object0$repeated, object1$repeated)) {
      stop("Repeated indices differs in the models", call. = FALSE)
    }
  }
  nms0 <- names(object0$coefficients)
  nms1 <- names(object1$coefficients)
  if (is.null(nms0) || is.null(nms1)) {
    stop("Models must have named coefficients to assess nesting", call. = FALSE)
  }
  if (length(nms0) == length(nms1)) {
    stop("Models must be nested (different numbers of coefficients expected)", call. = FALSE)
  }
  if (length(nms0) < length(nms1)) {
    obj0 <- object0
    obj1 <- object1
    small <- nms0
    big <- nms1
  } else {
    obj0 <- object1
    obj1 <- object0
    small <- nms1
    big <- nms0
  }
  if (length(setdiff(small, big)) != 0L) {
    stop("Models must be nested", call. = FALSE)
  }
  names_test <- setdiff(big, small)
  if (length(names_test) == 0L) {
    stop("Models must be nested", call. = FALSE)
  }
  index <- match(names_test, big)
  if (anyNA(index)) {
    stop("Internal error while matching coefficients for nesting test", call. = FALSE)
  }
  list(obj0 = obj0, obj1 = obj1, index = index)
}


## chisq-approximations for calculating the p-value for sum of independent
## chi sq random variables
lcsumchisq <- function(x, test_stat, pmethod = c("rao-scott", "satterthwaite")) {
  pmethod <- match.arg(pmethod)
  x <- Re(x)
  if (!is.numeric(test_stat) || length(test_stat) != 1L || !is.finite(test_stat)) {
    stop("'test_stat' must be a single finite numeric value", call. = FALSE)
  }
  if (!is.numeric(x) || length(x) < 1L || any(!is.finite(x))) {
    stop("'x' must be a non-empty numeric vector of finite values", call. = FALSE)
  }
  x_bar <- mean(x)
  test_df <- length(x)
  if (!is.finite(x_bar) || x_bar <= 0) {
    stop("invalid eigenvalues in chi-square mixture approximation", call. = FALSE)
  }
  if (pmethod == "rao-scott") {
    test_stat <- test_stat / x_bar
    test_p <- 1 - pchisq(test_stat, df = test_df)
  } else {
    alpha <- sum((x - x_bar)^2) / (test_df * x_bar^2)
    test_df <- test_df / (1 + alpha^2)
    test_stat <- test_stat / ((1 + alpha^2) * x_bar)
    test_p <- 1 - pchisq(test_stat, df = test_df)
  }

  list(test_stat = test_stat, test_df = test_df, test_p = test_p)
}


wald_test <- function(object0, object1,
                      cov_type = c("robust", "bias-corrected", "df-adjusted", "naive")) {
  cov_type <- match.arg(cov_type)
  nested_models <- check_nested_models(object0, object1)
  obj1 <- nested_models$obj1
  index <- nested_models$index
  test_df <- length(index)
  if (test_df < 1L) {
    stop("Wald test failed: no parameters to test (empty index set)", call. = FALSE)
  }
  coeffs_test <- as.numeric(obj1$coefficients[index])
  cov_mat <- vcov(obj1, cov_type = cov_type)
  cov_test <- cov_mat[index, index, drop = FALSE]
  test_stat <- tryCatch(
    as.numeric(crossprod(coeffs_test, solve(cov_test, coeffs_test))),
    error = function(e) stop("Wald test failed: covariance matrix is singular or invalid", call. = FALSE)
  )
  if (!is.finite(test_stat)) {
    stop("Wald test failed: non-finite test statistic", call. = FALSE)
  }
  if (test_stat < 0) {
    if (test_stat > -1e-10) {
      test_stat <- 0
    } else {
      stop("Wald test failed: negative test statistic", call. = FALSE)
    }
  }
  test_p <- 1 - pchisq(test_stat, df = test_df)
  list(test_stat = test_stat, test_df = test_df, test_p = test_p)
}


## working wald test
working_wald_test <- function(object0, object1,
                              cov_type = c("robust", "bias-corrected", "df-adjusted", "naive"),
                              pmethod = c("rao-scott", "satterthwaite")) {
  cov_type <- match.arg(cov_type)
  pmethod <- match.arg(pmethod)
  nested_models <- check_nested_models(object0, object1)
  obj1 <- nested_models$obj1
  index <- nested_models$index
  test_df <- length(index)
  if (test_df < 1L) {
    stop("Working Wald test failed: no parameters to test (empty index set)", call. = FALSE)
  }
  coeffs_test <- as.numeric(obj1$coefficients[index])
  naive_mat <- vcov(obj1, cov_type = "naive")
  cov_test <- naive_mat[index, index, drop = FALSE]
  stat_test <- tryCatch(
    as.numeric(crossprod(coeffs_test, solve(cov_test, coeffs_test))),
    error = function(e) stop("Working Wald test failed: naive covariance matrix is singular or invalid",
                             call. = FALSE)
  )
  if (!is.finite(stat_test)) {
    stop("Working Wald test failed: non-finite test statistic", call. = FALSE)
  }
  if (stat_test < 0) {
    if (stat_test > -1e-10) {
      stat_test <- 0
    } else {
      stop("Working Wald test failed: negative test statistic", call. = FALSE)
    }
  }
  robust_mat <- vcov(obj1, cov_type = cov_type)
  rob_test <- robust_mat[index, index, drop = FALSE]
  eigen_test <- tryCatch(
    eigen(solve(cov_test, rob_test), only.values = TRUE)$values,
    error = function(e) stop("Working Wald test failed: could not compute eigenvalues", call. = FALSE)
  )
  eigen_test <- Re(eigen_test)
  eigen_test[eigen_test < 0 & eigen_test > -1e-10] <- 0
  lcsumchisq(eigen_test, stat_test, pmethod)
}


## working likelihood ratio test
working_lr_test <- function(object0, object1,
                            cov_type = c("robust", "bias-corrected", "df-adjusted", "naive"),
                            pmethod = c("rao-scott", "satterthwaite")) {
  cov_type <- match.arg(cov_type)
  pmethod <- match.arg(pmethod)
  nested_models <- check_nested_models(object0, object1)
  obj1 <- nested_models$obj1
  obj0 <- nested_models$obj0
  index <- nested_models$index
  test_df <- length(index)
  if (test_df < 1L) {
    stop("Working LR test failed: no parameters to test (empty index set)", call. = FALSE)
  }
  if(obj0$phi != obj0$phi) {
    stop("Working LR test failed: dispersion parameters differ", call.=FALSE)
  }
  if(obj1$family %in% c("poisson", "binomial")) {
    if(obj0$phi != 1 || obj1$phi != 1)
      stop("Working LR test failed: dispersion parameters must be 1", call.=FALSE)
  }
  ll_obj0 <- compute_quasi_loglikelihood(obj0)
  ll_obj1 <- compute_quasi_loglikelihood(obj1)
  stat_test <- -2 * (ll_obj0 - ll_obj1)
  if (!is.finite(stat_test)) {
    stop("Working LR test failed: non-finite test statistic", call. = FALSE)
  }
  if (stat_test < 0) {
    if (stat_test > -1e-10) {
      stat_test <- 0
    } else {
      stop("Working LR test failed: negative test statistic", call. = FALSE)
    }
  }
  naive_mat <- vcov(obj1, cov_type = "naive")
  cov_test <- naive_mat[index, index, drop = FALSE]
  robust_mat <- vcov(obj1, cov_type = cov_type)
  rob_test <- robust_mat[index, index, drop = FALSE]
  eigen_test <- tryCatch(
    eigen(solve(cov_test, rob_test), only.values = TRUE)$values,
    error = function(e) stop("Working LR test failed: could not compute eigenvalues", call. = FALSE)
  )
  eigen_test <- Re(eigen_test)
  eigen_test[eigen_test < 0 & eigen_test > -1e-10] <- 0
  lcsumchisq(eigen_test, stat_test, pmethod)
}



## score test
score_test <- function(object0, object1,
                       cov_type = c("robust", "bias-corrected", "df-adjusted", "naive")) {
  cov_type <- match.arg(cov_type)
  nested_models <- check_nested_models(object0, object1)
  obj1 <- nested_models$obj1
  obj0 <- nested_models$obj0
  index <- nested_models$index
  test_df <- length(index)
  if (test_df < 1L) {
    stop("Score test failed: no parameters to test (empty index set)", call. = FALSE)
  }
  coeffs_test <- obj1$coefficients
  coeffs_test[index] <- 0
  coeffs_test[-index] <- obj0$coefficients
  is_geewa <- identical(as.character(obj1$call[[1]]), "geewa")
  if (is_geewa) {
    uvector <- estimating_equations_gee_cc(
      obj1$y, obj1$x, obj1$id, obj1$repeated, obj1$prior.weights,
      obj1$family$link, obj1$family$family,
      coeffs_test, obj0$fitted.values, obj0$linear.predictors,
      obj1$association_structure, obj1$alpha, obj1$phi
    )
    covariance <- get_covariance_matrices_cc(
      obj1$y, obj1$x, obj1$id, obj1$repeated, obj1$prior.weights,
      obj1$family$link, obj1$family$family,
      obj0$fitted.values, obj0$linear.predictors,
      obj1$association_structure, obj1$alpha, obj1$phi
    )
  } else {
    if (length(obj1$alpha) == 1) {
      association_alpha <- rep(obj1$alpha, choose(max(obj1$repeated), 2))
    } else {
      association_alpha <- obj1$alpha
    }
    uvector <- estimating_equations_gee_or(
      obj1$y, obj1$x, obj1$id, obj1$repeated, obj1$prior.weights,
      obj1$family$link,
      coeffs_test, obj0$fitted.values, obj0$linear.predictors,
      association_alpha
    )
    covariance <- get_covariance_matrices_or(
      obj1$y, obj1$x, obj1$id, obj1$repeated, obj1$prior.weights,
      obj1$family$link,
      obj0$fitted.values, obj0$linear.predictors,
      association_alpha
    )
  }
  cov_test <- switch(
    cov_type,
    robust = covariance$robust_covariance[index, index, drop = FALSE],
    naive = covariance$naive_covariance[index, index, drop = FALSE],
    `bias-corrected` = covariance$bc_covariance[index, index, drop = FALSE],
    {
      cl_no <- obj1$clusters_no
      coef_no <- length(coef(obj1))
      if (cl_no <= coef_no) {
        stop("Score test failed: clusters_no must be > number of coefficients for df-adjusted covariance",
             call. = FALSE)
      }
      (cl_no / (cl_no - coef_no)) * covariance$robust_covariance[index, index, drop = FALSE]
    }
  )
  naive_cov <- covariance$naive_covariance
  mid <- naive_cov[index, , drop = FALSE] %*% uvector
  sol <- tryCatch(
    solve(cov_test, mid),
    error = function(e) stop("Score test failed: covariance matrix is singular or invalid", call. = FALSE)
  )
  test_stat <- as.numeric((t(uvector) %*% naive_cov[, index, drop = FALSE]) %*% sol)
  if (!is.finite(test_stat)) {
    stop("Score test failed: non-finite test statistic", call. = FALSE)
  }
  if (test_stat < 0) {
    if (test_stat > -1e-10) {
      test_stat <- 0
    } else {
      stop("Score test failed: negative test statistic", call. = FALSE)
    }
  }
  test_p <- 1 - stats::pchisq(test_stat, df = test_df)
  list(test_stat = test_stat, test_df = test_df, test_p = test_p)
}



## working score test
working_score_test <- function(object0, object1,
                               cov_type = c("robust", "bias-corrected", "df-adjusted", "naive"),
                               pmethod = c("rao-scott", "satterthwaite")) {
  cov_type <- match.arg(cov_type)
  pmethod <- match.arg(pmethod)

  nested_models <- check_nested_models(object0, object1)
  obj1 <- nested_models$obj1
  obj0 <- nested_models$obj0
  index <- nested_models$index
  test_df <- length(index)

  if (test_df < 1L) {
    stop("Working score test failed: no parameters to test (empty index set)", call. = FALSE)
  }
  coeffs_test <- as.numeric(obj1$coefficients)
  names(coeffs_test) <- names(obj1$coefficients)
  coeffs_test[index] <- 0
  coeffs_test[-index] <- obj0$coefficients
  is_geewa <- identical(as.character(obj1$call[[1]]), "geewa")
  if (is_geewa) {
    uvector <- estimating_equations_gee_cc(
      obj1$y, obj1$x, obj1$id, obj1$repeated, obj1$prior.weights,
      obj1$family$link, obj1$family$family,
      coeffs_test, obj0$fitted.values, obj0$linear.predictors,
      obj1$association_structure, obj1$alpha, obj1$phi
    )
    covariance <- get_covariance_matrices_cc(
      obj1$y, obj1$x, obj1$id, obj1$repeated, obj1$prior.weights,
      obj1$family$link, obj1$family$family,
      obj0$fitted.values, obj0$linear.predictors,
      obj1$association_structure, obj1$alpha, obj1$phi
    )
  } else {
    if (length(obj1$alpha) == 1) {
      mrep <- max(as.integer(obj1$repeated))
      association_alpha <- rep(obj1$alpha, choose(mrep, 2))
    } else {
      association_alpha <- obj1$alpha
    }
    uvector <- estimating_equations_gee_or(
      obj1$y, obj1$x, obj1$id, obj1$repeated, obj1$prior.weights,
      obj1$family$link,
      coeffs_test, obj0$fitted.values, obj0$linear.predictors,
      association_alpha
    )

    covariance <- get_covariance_matrices_or(
      obj1$y, obj1$x, obj1$id, obj1$repeated, obj1$prior.weights,
      obj1$family$link,
      obj0$fitted.values, obj0$linear.predictors,
      association_alpha
    )
  }
  cov_test <- covariance$naive_covariance
  rob_test <- switch(
    cov_type,
    robust = covariance$robust_covariance,
    naive = covariance$naive_covariance,
    `bias-corrected` = covariance$bc_covariance,
    {
      m <- obj1$clusters_no
      p <- length(coef(obj1))
      if (!is.finite(m) || length(m) != 1L) {
        stop("Working score test failed: invalid clusters_no", call. = FALSE)
      }
      if (m <= p) {
        stop("Working score test failed: clusters_no must be > number of coefficients for df-adjusted covariance",
             call. = FALSE)
      }
      (m / (m - p)) * covariance$robust_covariance
    }
  )
  cov_ii <- cov_test[index, index, drop = FALSE]
  mid <- cov_test[index, , drop = FALSE] %*% uvector
  sol <- tryCatch(
    solve(cov_ii, mid),
    error = function(e) stop("Working score test failed: naive covariance submatrix is singular or invalid",
                             call. = FALSE)
  )
  score_stat <- as.numeric((t(uvector) %*% cov_test[, index, drop = FALSE]) %*% sol)
  if (!is.finite(score_stat)) {
    stop("Working score test failed: non-finite test statistic", call. = FALSE)
  }
  if (score_stat < 0) {
    if (score_stat > -1e-10) score_stat <- 0 else stop("Working score test failed: negative test statistic",
                                                       call. = FALSE)
  }
  eigen_test <- tryCatch(
    eigen(solve(cov_ii, rob_test[index, index, drop = FALSE]), only.values = TRUE)$values,
    error = function(e) stop("Working score test failed: could not compute eigenvalues", call. = FALSE)
  )
  eigen_test <- Re(eigen_test)
  eigen_test[eigen_test < 0 & eigen_test > -1e-10] <- 0
  lcsumchisq(eigen_test, score_stat, pmethod)
}



## anova list
anova_geerlist <- function(object, ..., test, cov_type, pmethod){
    response_vectors <-
      as.character(lapply(object, function(x) deparse(formula(x)[[2L]])))
    same_response <- response_vectors == response_vectors[1L]
    if (!all(same_response)) {
      object <- object[same_response]
      warning("Models with response ", deparse(response_vectors[!same_response]),
              " removed because response differs from Model 1")
    }
    if (test == "working-lrt") {
      association_vector <-
        as.character(lapply(object, function(x) x$association_structure))
      independence_vector <- association_vector == "independence"
      if (all(!independence_vector))
        stop("the modified working lr test requires independence working models")
      if (!all(independence_vector)) {
        object <- object[independence_vector]
        warning("Models with association structure ",
                deparse(association_vector[!independence_vector]),
                " removed because association structure differs from independence")
      }
    }
    association_vector <-
      as.character(lapply(object, function(x) x$association_structure))
    same_association <- association_vector == association_vector[1L]
    if (!all(same_association)) {
      object <- object[same_association]
      warning("Models with association structure ",
              deparse(association_vector[!same_association]),
              " removed because association structure differs from model 1")
    }
    sample_size_vector <- sapply(object, function(x) length(x$residuals))
    if (any(sample_size_vector != sample_size_vector[1L]))
      stop("models were not all fitted to the same size of dataset")
    models_no <- length(object)
    if (models_no == 1)
      return(anova.geer(object[[1L]], cov_type = cov_type, pmethod = pmethod))

    resdf  <- as.numeric(lapply(object, function(x) x$df.residual))
    table <- data.frame(resdf, c(NA, resdf[-1]), c(NA, resdf[-1]), c(NA, resdf[-1]))
    dimnames(table) <- list(1L:models_no,
                            c("Resid. Df", "Df", "Chi", "Pr(>Chi)"))
    test_type <- switch(test,
                        wald = "Wald",
                        score = "Score",
                        `working-wald` = "Working Wald",
                        `working-score` = "Working Score",
                        `working-lrt` = "Working LRT")
    for (i in 2:models_no) {
      value <-
        switch(test,
               wald =
                 wald_test(object[[i - 1]], object[[i]], cov_type),
               score =
                 score_test(object[[i - 1]], object[[i]], cov_type),
               `working-wald` =
                 working_wald_test(object[[i - 1]], object[[i]], cov_type, pmethod),
               `working-score` =
                 working_score_test(object[[i - 1]], object[[i]], cov_type, pmethod),
               `working-lrt` =
                 working_lr_test(object[[i - 1]], object[[i]], cov_type, pmethod)
        )
      table[i, -1] <- c(value$test_df, value$test_stat, value$test_p)
    }
    title  <- paste("Analysis of", test_type, "Statistic Table\n", sep = " ")
    variables <- lapply(object, function(x)
      paste(deparse(formula(x)), collapse = "\n") )
    topnote <-
      paste("Model ",
            format(1L:models_no),
            ": ",
            variables,
            sep = "",
            collapse = "\n")
    structure(table,
              heading = c(title, topnote),
              class = c("anova", "data.frame"))
  }
