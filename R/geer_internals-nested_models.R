## nested-model checks
check_nested_models <- function(object0, object1) {
  object0 <- check_geer_object(object0, "object0")
  object1 <- check_geer_object(object1, "object1")
  if (!is.null(object0$obs_no) && !is.null(object1$obs_no)) {
    if (!identical(object0$obs_no, object1$obs_no)) {
      stop("models were not fit on the same observations", call. = FALSE)
    }
  } else {
    if (!identical(object0$y, object1$y)) {
      stop("response variable differs between models", call. = FALSE)
    }
    if (!identical(object0$id, object1$id)) {
      stop("cluster identifiers differ between models", call. = FALSE)
    }
    if (!identical(object0$repeated, object1$repeated)) {
      stop("repeated indices differ between models", call. = FALSE)
    }
  }
  coef_names0 <- names(object0$coefficients)
  coef_names1 <- names(object1$coefficients)
  if (is.null(coef_names0) || is.null(coef_names1)) {
    stop("models must have named coefficients to assess nesting", call. = FALSE)
  }
  if (length(coef_names0) == length(coef_names1)) {
    stop("models must be nested and have different numbers of coefficients", call. = FALSE)
  }
  if (length(coef_names0) < length(coef_names1)) {
    smaller_model <- object0
    larger_model <- object1
    nms_small <- coef_names0
    nms_big <- coef_names1
  } else {
    smaller_model <- object1
    larger_model <- object0
    nms_small <- coef_names1
    nms_big <- coef_names0
  }
  if (length(setdiff(nms_small, nms_big)) != 0L) {
    stop("models must be nested", call. = FALSE)
  }
  index <- match(setdiff(nms_big, nms_small), nms_big)
  if (anyNA(index) || !length(index)) {
    stop("models must be nested", call. = FALSE)
  }
  list(
    obj0 = smaller_model,
    obj1 = larger_model,
    index = index
  )
}


## shared test helpers
.validate_test_index <- function(index, context) {
  test_df <- length(index)
  if (test_df < 1L) {
    stop(sprintf("%s failed: no parameters to test (empty index set)", context),
         call. = FALSE)
  }
  test_df
}


.validate_test_statistic <- function(test_stat, context, tol = 1e-10) {
  if (!is.finite(test_stat)) {
    stop(sprintf("%s failed: non-finite test statistic", context), call. = FALSE)
  }
  if (test_stat < 0) {
    if (test_stat > -tol) {
      test_stat <- 0
    } else {
      stop(sprintf("%s failed: negative test statistic", context), call. = FALSE)
    }
  }
  test_stat
}


.df_adjust_covariance <- function(robust_covariance,
                                  clusters_no,
                                  coef_no,
                                  context) {
  if (!is.finite(clusters_no) || length(clusters_no) != 1L) {
    stop(sprintf("%s failed: invalid clusters_no", context), call. = FALSE)
  }
  if (clusters_no <= coef_no) {
    stop(
      sprintf(
        "%s failed: clusters_no must be > number of coefficients for df-adjusted covariance",
        context
      ),
      call. = FALSE
    )
  }
  (clusters_no / (clusters_no - coef_no)) * robust_covariance
}


.compute_mixture_eigenvalues <- function(naive_mat, robust_mat, context) {
  eigenvalues <- tryCatch(
    eigen(solve(naive_mat, robust_mat), only.values = TRUE)$values,
    error = function(e) {
      stop(sprintf("%s failed: could not compute eigenvalues", context),
           call. = FALSE)
    }
  )
  eigenvalues <- Re(eigenvalues)
  eigenvalues[eigenvalues < 0 & eigenvalues > -1e-10] <- 0
  eigenvalues
}


.model_sample_size <- function(object) {
  if (!is.null(object$obs_no)) {
    return(as.integer(object$obs_no))
  }
  length(object$residuals)
}


.is_geewa_fit <- function(object) {
  identical(as.character(object$call[[1L]]), "geewa")
}


.get_or_alpha <- function(object) {
  if (length(object$alpha) == 1L) {
    rep(object$alpha, choose(max(as.integer(object$repeated)), 2))
  } else {
    object$alpha
  }
}


## chi-square mixture approximations
lcsumchisq <- function(x, test_stat, pmethod = c("rao-scott", "satterthwaite")) {
  pmethod <- match.arg(pmethod)
  x <- Re(x)
  if (!is.numeric(test_stat) || length(test_stat) != 1L || !is.finite(test_stat)) {
    stop("'test_stat' must be a single finite numeric value", call. = FALSE)
  }
  if (!is.numeric(x) || !length(x) || any(!is.finite(x))) {
    stop("'x' must be a non-empty numeric vector of finite values", call. = FALSE)
  }
  x_bar <- mean(x)
  test_df <- length(x)
  if (!is.finite(x_bar) || x_bar <= 0) {
    stop("invalid eigenvalues in chi-square mixture approximation", call. = FALSE)
  }
  if (identical(pmethod, "rao-scott")) {
    test_stat <- test_stat / x_bar
    test_p <- 1 - pchisq(test_stat, df = test_df)
  } else {
    satt_cv2 <- sum((x - x_bar)^2) / (test_df * x_bar^2)
    test_df <- test_df / (1 + satt_cv2)
    test_stat <- test_stat / ((1 + satt_cv2) * x_bar)
    test_p <- 1 - pchisq(test_stat, df = test_df)
  }

  list(test_stat = test_stat, test_df = test_df, test_p = test_p)
}


## score components
get_score_components <- function(obj0, obj1, test_coefficients) {
  if (.is_geewa_fit(obj1)) {
    score_vector <- estimating_equations_gee_cc(
      obj1$y, obj1$x, obj1$id, obj1$repeated, obj1$prior.weights,
      obj1$family$link, obj1$family$family,
      test_coefficients, obj0$fitted.values, obj0$linear.predictors,
      obj1$association_structure, obj1$alpha, obj1$phi
    )

    covariance <- get_covariance_matrices_cc(
      obj1$y, obj1$x, obj1$id, obj1$repeated, obj1$prior.weights,
      obj1$family$link, obj1$family$family,
      obj0$fitted.values, obj0$linear.predictors,
      obj1$association_structure, obj1$alpha, obj1$phi
    )
  } else {
    association_alpha <- .get_or_alpha(obj1)

    score_vector <- estimating_equations_gee_or(
      obj1$y, obj1$x, obj1$id, obj1$repeated, obj1$prior.weights,
      obj1$family$link,
      test_coefficients, obj0$fitted.values, obj0$linear.predictors,
      association_alpha
    )
    covariance <- get_covariance_matrices_or(
      obj1$y, obj1$x, obj1$id, obj1$repeated, obj1$prior.weights,
      obj1$family$link,
      obj0$fitted.values, obj0$linear.predictors,
      association_alpha
    )
  }
  list(
    score_vector = score_vector,
    naive_covariance = covariance$naive_covariance,
    robust_covariance = covariance$robust_covariance,
    bc_covariance = covariance$bc_covariance
  )
}


## Wald test
wald_test <- function(object0, object1,
                      cov_type = c("robust", "bias-corrected", "df-adjusted", "naive")) {
  cov_type <- match.arg(cov_type)
  nested_models <- check_nested_models(object0, object1)
  obj1 <- nested_models$obj1
  index <- nested_models$index
  test_df <- .validate_test_index(index, "Wald test")
  test_coefficients <- as.numeric(obj1$coefficients[index])
  cov_mat <- vcov(obj1, cov_type = cov_type)
  cov_test <- cov_mat[index, index, drop = FALSE]
  test_stat <- tryCatch(
    as.numeric(crossprod(test_coefficients, solve(cov_test, test_coefficients))),
    error = function(e) {
      stop("Wald test failed: covariance matrix is singular or invalid",
           call. = FALSE)
    }
  )
  test_stat <- .validate_test_statistic(test_stat, "Wald test")
  test_p <- 1 - pchisq(test_stat, df = test_df)
  list(test_stat = test_stat, test_df = test_df, test_p = test_p)
}


## working Wald test
working_wald_test <- function(object0, object1,
                              cov_type = c("robust", "bias-corrected", "df-adjusted", "naive"),
                              pmethod = c("rao-scott", "satterthwaite")) {
  cov_type <- match.arg(cov_type)
  pmethod <- match.arg(pmethod)
  nested_models <- check_nested_models(object0, object1)
  obj1 <- nested_models$obj1
  index <- nested_models$index
  .validate_test_index(index, "Working Wald test")
  test_coefficients <- as.numeric(obj1$coefficients[index])
  naive_mat <- vcov(obj1, cov_type = "naive")
  cov_test <- naive_mat[index, index, drop = FALSE]
  stat_test <- tryCatch(
    as.numeric(crossprod(test_coefficients, solve(cov_test, test_coefficients))),
    error = function(e) {
      stop(
        "Working Wald test failed: naive covariance matrix is singular or invalid",
        call. = FALSE
      )
    }
  )
  stat_test <- .validate_test_statistic(stat_test, "Working Wald test")
  robust_mat <- vcov(obj1, cov_type = cov_type)
  rob_test <- robust_mat[index, index, drop = FALSE]
  eigen_test <- .compute_mixture_eigenvalues(
    naive_mat = cov_test,
    robust_mat = rob_test,
    context = "Working Wald test"
  )
  lcsumchisq(eigen_test, stat_test, pmethod)
}


## working likelihood ratio test
working_lr_test <- function(object0, object1,
                            cov_type = c("robust", "bias-corrected", "df-adjusted", "naive"),
                            pmethod = c("rao-scott", "satterthwaite")) {
  cov_type <- match.arg(cov_type)
  pmethod <- match.arg(pmethod)
  nested_models <- check_nested_models(object0, object1)
  obj0 <- nested_models$obj0
  obj1 <- nested_models$obj1
  index <- nested_models$index
  .validate_test_index(index, "Working LR test")
  phi_tol <- sqrt(.Machine$double.eps) * max(abs(obj0$phi), abs(obj1$phi), 1)
  if (abs(obj0$phi - obj1$phi) > phi_tol) {
    stop("Working LR test failed: dispersion parameters differ between models",
         call. = FALSE)
  }
  if (obj1$family$family %in% c("poisson", "binomial")) {
    if (abs(obj0$phi - 1) > phi_tol || abs(obj1$phi - 1) > phi_tol) {
      stop(
        "Working LR test failed: dispersion parameter must equal 1 for Poisson/binomial models",
        call. = FALSE
      )
    }
  }
  ll_obj0 <- compute_quasi_loglikelihood(obj0)
  ll_obj1 <- compute_quasi_loglikelihood(obj1)
  stat_test <- -2 * (ll_obj0 - ll_obj1)
  stat_test <- .validate_test_statistic(stat_test, "Working LR test")
  naive_mat <- vcov(obj1, cov_type = "naive")
  cov_test <- naive_mat[index, index, drop = FALSE]
  robust_mat <- vcov(obj1, cov_type = cov_type)
  rob_test <- robust_mat[index, index, drop = FALSE]
  eigen_test <- .compute_mixture_eigenvalues(
    naive_mat = cov_test,
    robust_mat = rob_test,
    context = "Working LR test"
  )
  lcsumchisq(eigen_test, stat_test, pmethod)
}


## score test
score_test <- function(object0, object1,
                       cov_type = c("robust", "bias-corrected", "df-adjusted", "naive")) {
  cov_type <- match.arg(cov_type)
  nested_models <- check_nested_models(object0, object1)
  obj0 <- nested_models$obj0
  obj1 <- nested_models$obj1
  index <- nested_models$index
  test_df <- .validate_test_index(index, "Score test")
  test_coefficients <- obj1$coefficients
  test_coefficients[index] <- 0
  test_coefficients[names(obj0$coefficients)] <- obj0$coefficients
  sc <- get_score_components(obj0, obj1, test_coefficients)
  score_vector <- sc$score_vector
  cov_test <- switch(
    cov_type,
    robust = sc$robust_covariance[index, index, drop = FALSE],
    naive = sc$naive_covariance[index, index, drop = FALSE],
    `bias-corrected` = sc$bc_covariance[index, index, drop = FALSE],
    `df-adjusted` = .df_adjust_covariance(
      robust_covariance = sc$robust_covariance[index, index, drop = FALSE],
      clusters_no = obj1$clusters_no,
      coef_no = length(coef(obj1)),
      context = "Score test"
    )
  )
  naive_cov <- sc$naive_covariance
  mid <- naive_cov[index, , drop = FALSE] %*% score_vector
  sol <- tryCatch(
    solve(cov_test, mid),
    error = function(e) {
      stop("Score test failed: covariance matrix is singular or invalid",
           call. = FALSE)
    }
  )
  test_stat <- as.numeric((t(score_vector) %*% naive_cov[, index, drop = FALSE]) %*% sol)
  test_stat <- .validate_test_statistic(test_stat, "Score test")
  test_p <- 1 - pchisq(test_stat, df = test_df)
  list(test_stat = test_stat, test_df = test_df, test_p = test_p)
}


## working score test -------------------------------------------------------
working_score_test <- function(object0, object1,
                               cov_type = c("robust", "bias-corrected", "df-adjusted", "naive"),
                               pmethod = c("rao-scott", "satterthwaite")) {
  cov_type <- match.arg(cov_type)
  pmethod <- match.arg(pmethod)
  nested_models <- check_nested_models(object0, object1)
  obj0 <- nested_models$obj0
  obj1 <- nested_models$obj1
  index <- nested_models$index
  .validate_test_index(index, "Working score test")
  test_coefficients <- as.numeric(obj1$coefficients)
  names(test_coefficients) <- names(obj1$coefficients)
  test_coefficients[index] <- 0
  test_coefficients[names(obj0$coefficients)] <- obj0$coefficients
  sc <- get_score_components(obj0, obj1, test_coefficients)
  score_vector <- sc$score_vector
  cov_test <- sc$naive_covariance
  rob_test <- switch(
    cov_type,
    robust = sc$robust_covariance,
    naive = sc$naive_covariance,
    `bias-corrected` = sc$bc_covariance,
    `df-adjusted` = .df_adjust_covariance(
      robust_covariance = sc$robust_covariance,
      clusters_no = obj1$clusters_no,
      coef_no = length(coef(obj1)),
      context = "Working score test"
    )
  )
  cov_ii <- cov_test[index, index, drop = FALSE]
  mid <- cov_test[index, , drop = FALSE] %*% score_vector
  sol <- tryCatch(
    solve(cov_ii, mid),
    error = function(e) {
      stop(
        "Working score test failed: naive covariance submatrix is singular or invalid",
        call. = FALSE
      )
    }
  )
  score_stat <- as.numeric((t(score_vector) %*% cov_test[, index, drop = FALSE]) %*% sol)
  score_stat <- .validate_test_statistic(score_stat, "Working score test")
  eigen_test <- .compute_mixture_eigenvalues(
    naive_mat = cov_ii,
    robust_mat = rob_test[index, index, drop = FALSE],
    context = "Working score test"
  )
  lcsumchisq(eigen_test, score_stat, pmethod)
}


## anova for a list of geer fits
anova_geerlist <- function(object, ..., test, cov_type, pmethod) {
  response_vector <- vapply(
    object,
    function(x) paste(deparse(formula(x)[[2L]]), collapse = " "),
    character(1)
  )
  same_response <- response_vector == response_vector[1L]
  if (!all(same_response)) {
    removed <- unique(response_vector[!same_response])
    object <- object[same_response]
    warning(
      "models with response ",
      paste(removed, collapse = ", "),
      " removed because response differs from Model 1",
      call. = FALSE
    )
  }
  if (identical(test, "working-lrt")) {
    association_vector <- vapply(
      object,
      function(x) x$association_structure,
      character(1)
    )
    independence_vector <- association_vector == "independence"

    if (all(!independence_vector)) {
      stop("the modified working LRT requires independence working models",
           call. = FALSE)
    }

    if (!all(independence_vector)) {
      removed <- unique(association_vector[!independence_vector])
      object <- object[independence_vector]
      warning(
        "models with association structure ",
        paste(removed, collapse = ", "),
        " removed because association structure differs from independence",
        call. = FALSE
      )
    }
  }
  association_vector <- vapply(
    object,
    function(x) x$association_structure,
    character(1)
  )
  same_association <- association_vector == association_vector[1L]
  if (!all(same_association)) {
    removed <- unique(association_vector[!same_association])
    object <- object[same_association]
    warning(
      "models with association structure ",
      paste(removed, collapse = ", "),
      " removed because association structure differs from Model 1",
      call. = FALSE
    )
  }
  sample_size_vector <- vapply(object, .model_sample_size, integer(1))
  if (any(sample_size_vector != sample_size_vector[1L])) {
    stop("models were not all fitted to the same size of dataset", call. = FALSE)
  }
  models_no <- length(object)
  if (models_no == 1L) {
    return(
      anova.geer(
        object[[1L]],
        test = test,
        cov_type = cov_type,
        pmethod = pmethod
      )
    )
  }
  resdf <- vapply(object, function(x) x$df.residual, numeric(1))
  table <- data.frame(
    `Resid. Df` = resdf,
    Df = NA_real_,
    Chi = NA_real_,
    `Pr(>Chi)` = NA_real_,
    check.names = FALSE
  )
  test_type <- test_label(test)
  for (i in 2:models_no) {
    value <- switch(
      test,
      wald = wald_test(object[[i - 1L]], object[[i]], cov_type),
      score = score_test(object[[i - 1L]], object[[i]], cov_type),
      `working-wald` = working_wald_test(object[[i - 1L]], object[[i]], cov_type, pmethod),
      `working-score` = working_score_test(object[[i - 1L]], object[[i]], cov_type, pmethod),
      `working-lrt` = working_lr_test(object[[i - 1L]], object[[i]], cov_type, pmethod)
    )
    table[i, -1L] <- c(value$test_df, value$test_stat, value$test_p)
  }
  title <- paste("Analysis of", test_type, "Statistic Table\n")
  variables <- vapply(
    object,
    function(x) paste(deparse(formula(x)), collapse = "\n"),
    character(1)
  )
  topnote <- paste(
    "Model ",
    format(seq_len(models_no)),
    ": ",
    variables,
    sep = "",
    collapse = "\n"
  )
  structure(
    table,
    heading = c(title, topnote),
    class = c("anova", "data.frame")
  )
}
