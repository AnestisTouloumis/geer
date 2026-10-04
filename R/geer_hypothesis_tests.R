## nested-model checks
## Settings that must agree for two fits to be compared. Nesting by coefficient
## names alone does not guarantee that the two models differ only in their mean
## structure: a Wald or score comparison of, say, a bcgee-robust fit with a
## brgee-naive fit, or of fits with different dispersion handling, would
## otherwise run silently.
check_comparable_fit_settings <- function(object0, object1) {
  same_or_stop <- function(value0, value1, what) {
    if (!identical(value0, value1)) {
      stop("models differ in ", what, call. = FALSE)
    }
  }
  same_or_stop(get_geer_fit_function(object0), get_geer_fit_function(object1),
               "the fitting function ('geewa' or 'geewa_binary')")
  same_or_stop(object0$method, object1$method, "the estimation method")
  same_or_stop(object0$association_structure, object1$association_structure,
               "the working association structure")
  structure_name <- object0$association_structure
  if (identical(structure_name, "m-dependent")) {
    same_or_stop(length(object0$alpha), length(object1$alpha),
                 "the order of the m-dependent working structure")
  }
  if (identical(structure_name, "fixed")) {
    same_or_stop(as.numeric(object0$alpha), as.numeric(object1$alpha),
                 "the fixed association parameters")
  }
  if (identical(get_geer_fit_function(object0), "geewa")) {
    phi_fixed0 <- is_phi_fixed(object0)
    phi_fixed1 <- is_phi_fixed(object1)
    same_or_stop(phi_fixed0, phi_fixed1, "the treatment of the dispersion parameter")
    if (phi_fixed0) {
      same_or_stop(as.numeric(object0$phi), as.numeric(object1$phi),
                   "the fixed dispersion value")
    }
    use_p0 <- if (is.null(object0$use_p)) TRUE else isTRUE(object0$use_p)
    use_p1 <- if (is.null(object1$use_p)) TRUE else isTRUE(object1$use_p)
    same_or_stop(use_p0, use_p1, "'use_p'")
  }
  contrasts0 <- object0$contrasts
  contrasts1 <- object1$contrasts
  shared <- intersect(names(contrasts0), names(contrasts1))
  if (length(shared)) {
    same_or_stop(contrasts0[shared], contrasts1[shared], "the contrasts")
  }
  invisible(TRUE)
}


check_nested_models <- function(object0, object1) {
  object0 <- check_geer_object(object0, "object0")
  object1 <- check_geer_object(object1, "object1")
  if (!identical(object0$obs_no, object1$obs_no)) {
    stop("models were not fit on the same number of observations", call. = FALSE)
  }
  if (!identical(object0$y, object1$y)) {
    stop("response variable differs between models", call. = FALSE)
  }
  if (!identical(object0$id, object1$id)) {
    stop("cluster identifiers differ between models", call. = FALSE)
  }
  if (!identical(object0$repeated, object1$repeated)) {
    stop("repeated indices differ between models", call. = FALSE)
  }
  if (!identical(object0$prior.weights, object1$prior.weights)) {
    stop("prior weights differ between models", call. = FALSE)
  }
  if (!identical(object0$offset, object1$offset)) {
    stop("offset differs between models", call. = FALSE)
  }
  if (!identical(object0$family$family, object1$family$family) ||
      !identical(object0$family$link, object1$family$link)) {
    stop("family or link differs between models", call. = FALSE)
  }
  check_comparable_fit_settings(object0, object1)
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
    object0 = smaller_model,
    object1 = larger_model,
    index = index
  )
}


check_test_index <- function(index, context) {
  test_df <- length(index)
  if (test_df < 1L) {
    stop(sprintf("%s failed: no parameters to test (empty index set)", context),
         call. = FALSE)
  }
  test_df
}


check_test_statistic <- function(test_stat, context, tol = 1e-10) {
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


compute_df_adjusted_covariance <- function(robust_covariance,
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


compute_mixture_eigenvalues <- function(naive_mat, robust_mat, context) {
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


get_model_sample_size <- function(object) {
  if (!is.null(object$obs_no)) {
    return(as.integer(object$obs_no))
  }
  length(object$residuals)
}


## The fitting function is recorded in the fit. The call is parsed only as a
## fallback for objects that do not carry it, and then only its last component,
## so that geer::geewa(...) is recognized. Parsing the call alone is not
## reliable: the head of geer::geewa(...) is `::`, and do.call(geewa, ...)
## stores the function object itself.
get_geer_fit_function <- function(object) {
  fit_function <- object$fit_function
  if (is.character(fit_function) && length(fit_function) == 1L &&
      !is.na(fit_function)) {
    return(fit_function)
  }
  call_head <- if (is.null(object$call)) NULL else object$call[[1L]]
  if (is.name(call_head) || is.call(call_head)) {
    call_head <- as.character(call_head)
    call_head <- call_head[[length(call_head)]]
    if (call_head %in% c("geewa", "geewa_binary")) {
      return(call_head)
    }
  }
  NA_character_
}


is_geewa_fit <- function(object) {
  fit_function <- get_geer_fit_function(object)
  if (is.na(fit_function)) {
    stop(
      "cannot determine whether the model was fitted by 'geewa' or 'geewa_binary'",
      call. = FALSE
    )
  }
  identical(fit_function, "geewa")
}


get_alpha_or <- function(object) {
  if (length(object$alpha) == 1L) {
    rep(object$alpha, choose(max(as.integer(object$repeated)), 2))
  } else {
    object$alpha
  }
}


compute_chisq_mixture <- function(x, test_stat, pmethod = geer_pmethod_choices) {
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
    test_p <- stats::pchisq(test_stat, df = test_df, lower.tail = FALSE)
  } else {
    satt_cv2 <- sum((x - x_bar)^2) / (test_df * x_bar^2)
    test_df <- test_df / (1 + satt_cv2)
    test_stat <- test_stat / ((1 + satt_cv2) * x_bar)
    test_p <- stats::pchisq(test_stat, df = test_df, lower.tail = FALSE)
  }
  list(test_stat = test_stat, test_df = test_df, test_p = test_p)
}


compute_score_components <- function(object0, object1) {
  if (is_geewa_fit(object1)) {
    score_vector <- estimating_equations_gee_cc(
      object1$y, object1$x, object1$id, object1$repeated, object1$prior.weights,
      object1$family$link, object1$family$family,
      object0$fitted.values, object0$linear.predictors,
      object1$association_structure, object1$alpha, object1$phi
    )
    covariance <- get_covariance_matrices_cc(
      object1$y, object1$x, object1$id, object1$repeated, object1$prior.weights,
      object1$family$link, object1$family$family,
      object0$fitted.values, object0$linear.predictors,
      object1$association_structure, object1$alpha, object1$phi
    )
  } else {
    association_alpha <- get_alpha_or(object1)
    score_vector <- estimating_equations_gee_or(
      object1$y, object1$x, object1$id, object1$repeated, object1$prior.weights,
      object1$family$link,
      object0$fitted.values, object0$linear.predictors,
      association_alpha
    )
    covariance <- get_covariance_matrices_or(
      object1$y, object1$x, object1$id, object1$repeated, object1$prior.weights,
      object1$family$link,
      object0$fitted.values, object0$linear.predictors,
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


compute_wald_test <- function(object0, object1,
                      cov_type = geer_cov_type_choices) {
  cov_type <- match.arg(cov_type)
  nested_models <- check_nested_models(object0, object1)
  obj1 <- nested_models$object1
  index <- nested_models$index
  test_df <- check_test_index(index, "Wald test")
  test_coefficients <- as.numeric(obj1$coefficients[index])
  cov_test <- stats::vcov(obj1, cov_type = cov_type)[index, index, drop = FALSE]
  test_stat <- tryCatch(
    as.numeric(crossprod(test_coefficients, solve(cov_test, test_coefficients))),
    error = function(e) {
      stop("Wald test failed: covariance matrix is singular or invalid",
           call. = FALSE)
    }
  )
  test_stat <- check_test_statistic(test_stat, "Wald test")
  test_p <- stats::pchisq(test_stat, df = test_df, lower.tail = FALSE)
  list(test_stat = test_stat, test_df = test_df, test_p = test_p)
}


compute_working_wald_test <- function(object0, object1,
                              cov_type = geer_cov_type_choices,
                              pmethod = geer_pmethod_choices) {
  cov_type <- match.arg(cov_type)
  pmethod <- match.arg(pmethod)
  nested_models <- check_nested_models(object0, object1)
  obj1 <- nested_models$object1
  index <- nested_models$index
  check_test_index(index, "Working Wald test")
  test_coefficients <- as.numeric(obj1$coefficients[index])
  naive_mat <- stats::vcov(obj1, cov_type = "naive")
  cov_test <- naive_mat[index, index, drop = FALSE]
  test_stat <- tryCatch(
    as.numeric(crossprod(test_coefficients, solve(cov_test, test_coefficients))),
    error = function(e) {
      stop(
        "Working Wald test failed: naive covariance matrix is singular or invalid",
        call. = FALSE
      )
    }
  )
  test_stat <- check_test_statistic(test_stat, "Working Wald test")
  robust_mat <- stats::vcov(obj1, cov_type = cov_type)
  rob_test <- robust_mat[index, index, drop = FALSE]
  eigen_test <- compute_mixture_eigenvalues(
    naive_mat = cov_test,
    robust_mat = rob_test,
    context = "Working Wald test"
  )
  compute_chisq_mixture(eigen_test, test_stat, pmethod)
}


compute_working_lrt_test <- function(object0, object1,
                             cov_type = geer_cov_type_choices,
                             pmethod = geer_pmethod_choices) {
  cov_type <- match.arg(cov_type)
  pmethod <- match.arg(pmethod)
  nested_models <- check_nested_models(object0, object1)
  obj0 <- nested_models$object0
  obj1 <- nested_models$object1
  index <- nested_models$index
  check_test_index(index, "Working LR test")
  if (obj1$family$family %in% c("poisson", "binomial")) {
    phi_tol <- 1e-6
    if (abs(obj0$phi - 1) > phi_tol || abs(obj1$phi - 1) > phi_tol) {
      stop(
        "Working LR test failed: dispersion parameter must equal 1 for Poisson/binomial models",
        call. = FALSE
      )
    }
  }
  ll_obj0 <- compute_quasi_loglikelihood(obj0)
  ll_obj1 <- compute_quasi_loglikelihood(obj1)
  test_stat <- -2 * (ll_obj0 - ll_obj1)
  test_stat <- check_test_statistic(test_stat, "Working LR test")
  naive_mat <- stats::vcov(obj1, cov_type = "naive")
  cov_test <- naive_mat[index, index, drop = FALSE]
  robust_mat <- stats::vcov(obj1, cov_type = cov_type)
  rob_test <- robust_mat[index, index, drop = FALSE]
  eigen_test <- compute_mixture_eigenvalues(
    naive_mat = cov_test,
    robust_mat = rob_test,
    context = "Working LR test"
  )
  compute_chisq_mixture(eigen_test, test_stat, pmethod)
}


compute_score_test <- function(object0, object1,
                       cov_type = geer_cov_type_choices) {
  cov_type <- match.arg(cov_type)
  nested_models <- check_nested_models(object0, object1)
  obj0 <- nested_models$object0
  obj1 <- nested_models$object1
  index <- nested_models$index
  test_df <- check_test_index(index, "Score test")
  sc <- compute_score_components(obj0, obj1)
  score_vector <- sc$score_vector
  cov_test <- switch(
    cov_type,
    robust = sc$robust_covariance[index, index, drop = FALSE],
    naive = sc$naive_covariance[index, index, drop = FALSE],
    `bias-corrected` = sc$bc_covariance[index, index, drop = FALSE],
    `df-adjusted` = compute_df_adjusted_covariance(
      robust_covariance = sc$robust_covariance[index, index, drop = FALSE],
      clusters_no = obj1$clusters_no,
      coef_no = length(stats::coef(obj1)),
      context = "Score test"
    ),
    jackknife = stats::vcov(obj1, cov_type = "jackknife")[index, index, drop = FALSE]
  )
  naive_mat <- sc$naive_covariance
  mid <- naive_mat[index, , drop = FALSE] %*% score_vector
  sol <- tryCatch(
    solve(cov_test, mid),
    error = function(e) {
      stop("Score test failed: covariance matrix is singular or invalid",
           call. = FALSE)
    }
  )
  test_stat <- as.numeric((t(score_vector) %*% naive_mat[, index, drop = FALSE]) %*% sol)
  test_stat <- check_test_statistic(test_stat, "Score test")
  test_p <- stats::pchisq(test_stat, df = test_df, lower.tail = FALSE)
  list(test_stat = test_stat, test_df = test_df, test_p = test_p)
}


compute_working_score_test <- function(object0, object1,
                               cov_type = geer_cov_type_choices,
                               pmethod = geer_pmethod_choices) {
  cov_type <- match.arg(cov_type)
  pmethod <- match.arg(pmethod)
  nested_models <- check_nested_models(object0, object1)
  obj0 <- nested_models$object0
  obj1 <- nested_models$object1
  index <- nested_models$index
  check_test_index(index, "Working score test")
  sc <- compute_score_components(obj0, obj1)
  score_vector <- sc$score_vector
  naive_mat <- sc$naive_covariance
  robust_mat <- switch(
    cov_type,
    robust = sc$robust_covariance,
    naive = sc$naive_covariance,
    `bias-corrected` = sc$bc_covariance,
    `df-adjusted` = compute_df_adjusted_covariance(
      robust_covariance = sc$robust_covariance,
      clusters_no = obj1$clusters_no,
      coef_no = length(stats::coef(obj1)),
      context = "Working score test"
    ),
    jackknife = stats::vcov(obj1, cov_type = "jackknife")
  )
  cov_test <- naive_mat[index, index, drop = FALSE]
  mid <- naive_mat[index, , drop = FALSE] %*% score_vector
  sol <- tryCatch(
    solve(cov_test, mid),
    error = function(e) {
      stop(
        "Working score test failed: naive covariance submatrix is singular or invalid",
        call. = FALSE
      )
    }
  )
  test_stat <- as.numeric((t(score_vector) %*% naive_mat[, index, drop = FALSE]) %*% sol)
  test_stat <- check_test_statistic(test_stat, "Working score test")
  eigen_test <- compute_mixture_eigenvalues(
    naive_mat = cov_test,
    robust_mat = robust_mat[index, index, drop = FALSE],
    context = "Working score test"
  )
  compute_chisq_mixture(eigen_test, test_stat, pmethod)
}


compute_anova_geer_list <- function(object, ..., test, cov_type, pmethod) {
  response_vector <- vapply(
    object,
    function(x) paste(deparse(stats::formula(x)[[2L]]), collapse = " "),
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
    if (!all(independence_vector)) {
      stop(
        "the modified working LRT requires all models to use an independence working structure",
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
  sample_size_vector <- vapply(object, get_model_sample_size, integer(1))
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
  test_type <- format_test_label(test)
  for (i in 2:models_no) {
    value <- switch(
      test,
      wald = compute_wald_test(object[[i - 1L]], object[[i]], cov_type),
      score = compute_score_test(object[[i - 1L]], object[[i]], cov_type),
      `working-wald` = compute_working_wald_test(object[[i - 1L]], object[[i]], cov_type, pmethod),
      `working-score` = compute_working_score_test(object[[i - 1L]], object[[i]], cov_type, pmethod),
      `working-lrt` = compute_working_lrt_test(object[[i - 1L]], object[[i]], cov_type, pmethod)
    )
    table[i, -1L] <- c(value$test_df, value$test_stat, value$test_p)
  }
  title <- paste("Analysis of", test_type, "Statistic Table\n")
  variables <- vapply(
    object,
    function(x) paste(deparse(stats::formula(x)), collapse = "\n"),
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
    class = c("geer_anova", "anova", "data.frame")
  )
}
