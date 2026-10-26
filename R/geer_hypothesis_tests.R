check_comparable_fit_settings <- function(object0, object1) {
  same_or_stop <- function(value0, value1, what) {
    if (!identical(value0, value1)) {
      stop("models differ in ", what, call. = FALSE)
    }
  }
  same_or_stop(
    get_geer_fit_function(object0),
    get_geer_fit_function(object1),
    "the fitting function (geewa() or geewa_binary())"
  )
  same_or_stop(object0$method, object1$method, "the estimation method")
  same_or_stop(
    object0$association_structure,
    object1$association_structure,
    "the working association structure"
  )
  structure_name <- object0$association_structure
  if (identical(structure_name, "m-dependent")) {
    same_or_stop(
      length(object0$alpha),
      length(object1$alpha),
      "the order of the m-dependent working structure"
    )
  }
  if (identical(structure_name, "fixed")) {
    same_or_stop(
      as.numeric(object0$alpha),
      as.numeric(object1$alpha),
      "the fixed association parameters"
    )
  }
  if (identical(get_geer_fit_function(object0), "geewa")) {
    phi_fixed0 <- is_phi_fixed(object0)
    phi_fixed1 <- is_phi_fixed(object1)
    same_or_stop(
      phi_fixed0,
      phi_fixed1,
      "the treatment of the dispersion parameter"
    )
    if (phi_fixed0) {
      same_or_stop(
        as.numeric(object0$phi),
        as.numeric(object1$phi),
        "the fixed dispersion value"
      )
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
    stop(
      "models must be fitted to the same observations",
      call. = FALSE
    )
  }
  if (!identical(object0$y, object1$y)) {
    stop("models differ in the response", call. = FALSE)
  }
  if (!identical(object0$id, object1$id)) {
    stop("models differ in the cluster identifiers", call. = FALSE)
  }
  if (!identical(object0$repeated, object1$repeated)) {
    stop("models differ in the repeated indices", call. = FALSE)
  }
  if (!identical(object0$prior.weights, object1$prior.weights)) {
    stop("models differ in the prior weights", call. = FALSE)
  }
  if (!identical(object0$offset, object1$offset)) {
    stop("models differ in the offset", call. = FALSE)
  }
  if (!identical(object0$family$family, object1$family$family) ||
      !identical(object0$family$link, object1$family$link)) {
    stop("models differ in the family or link", call. = FALSE)
  }
  check_comparable_fit_settings(object0, object1)
  coef_names0 <- names(object0$coefficients)
  coef_names1 <- names(object1$coefficients)
  if (is.null(coef_names0) || is.null(coef_names1)) {
    stop("models must have named coefficients to assess nesting", call. = FALSE)
  }
  if (length(coef_names0) == length(coef_names1)) {
    stop(
      "models must be nested and have different numbers of coefficients",
      call. = FALSE
    )
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
  check_nested_designs(smaller_model, larger_model)
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


check_nested_designs <- function(smaller_model, larger_model) {
  shared <- names(smaller_model$coefficients)
  design <- cbind(
    larger_model$x[, shared, drop = FALSE],
    smaller_model$x
  )
  scales <- vapply(
    seq_len(ncol(design)),
    function(index) max(abs(design[, index])),
    numeric(1)
  )
  if (any(!is.finite(design)) || any(scales == 0) ||
        qr(sweep(design, 2L, scales, `/`))$rank != length(shared)) {
    stop("models must be nested", call. = FALSE)
  }
  invisible(NULL)
}


check_test_index <- function(index, context) {
  test_df <- length(index)
  if (test_df < 1L) {
    stop(
      sprintf("%s failed: the models have no coefficients to test", context),
      call. = FALSE
    )
  }
  test_df
}


check_test_statistic <- function(test_stat, context, test_label,
                                 tol = 1e-10) {
  if (!is.finite(test_stat)) {
    stop(
      sprintf("%s failed: the %s statistic is not finite", context, test_label),
      call. = FALSE
    )
  }
  if (test_stat < 0) {
    if (test_stat <= -tol) {
      stop(
        sprintf(
          paste0(
            "%s failed: the %s statistic is negative (%s), probably because ",
            "the covariance estimate is not positive semi-definite"
          ),
          context, test_label, format(test_stat, digits = 3)
        ),
        call. = FALSE
      )
    }
    test_stat <- 0
  }
  test_stat
}


compute_df_adjusted_covariance <- function(robust_covariance,
                                           clusters_no,
                                           coef_no,
                                           context) {
  clusters_no <- check_integer_at_least(clusters_no, "clusters_no")
  if (clusters_no <= coef_no) {
    stop(
      sprintf(
        paste0(
          "%s failed: the df-adjusted covariance requires more clusters than ",
          "coefficients"
        ),
        context
      ),
      call. = FALSE
    )
  }
  (clusters_no / (clusters_no - coef_no)) * robust_covariance
}


check_bc_covariance <- function(clusters_no, coef_no, context) {
  if (clusters_no <= coef_no) {
    stop(
      sprintf(
        paste0(
          "%s failed: the bias-corrected covariance requires more clusters ",
          "than coefficients; use cov_type = \"robust\" or \"jackknife\""
        ),
        context
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}


compute_mixture_eigenvalues <- function(naive_mat, robust_mat, context,
                                        test_label) {
  eigenvalues <- tryCatch(
    eigen(solve(naive_mat, robust_mat), only.values = TRUE)$values,
    error = function(e) {
      stop(
        sprintf(
          paste0(
            "%s failed: the eigenvalues of the %s weight matrix could not be ",
            "computed"
          ),
          context, test_label
        ),
        call. = FALSE
      )
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
    stop("'object' must be a geewa() or geewa_binary() fit", call. = FALSE)
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


compute_chisq_mixture <- function(x,
                                  test_stat,
                                  pmethod,
                                  context) {
  x <- Re(x)
  if (!is.numeric(test_stat) || length(test_stat) != 1L ||
      !is.finite(test_stat)) {
    stop("'test_stat' must be a finite number", call. = FALSE)
  }
  if (!is.numeric(x) || !length(x) || any(!is.finite(x))) {
    stop(
      "'x' must be a non-empty numeric vector of finite values",
      call. = FALSE
    )
  }
  x_bar <- mean(x)
  test_df <- length(x)
  if (!is.finite(x_bar) || x_bar <= 0) {
    stop(
      sprintf(
        paste0(
          "%s failed: the weights of the chi-squared mixture are not finite ",
          "and positive"
        ),
        context
      ),
      call. = FALSE
    )
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
    score_vector <- get_estimating_equations_cc(
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
    score_vector <- get_estimating_equations_or(
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
    score_vector = as.numeric(score_vector),
    naive_covariance = covariance$naive_covariance,
    robust_covariance = covariance$robust_covariance,
    bc_covariance = covariance$bc_covariance
  )
}


compute_score_ss_cov <- function(object0,
                                 object1,
                                 cov_type,
                                 context) {
  select_small_sample_covariance(
    covariances = compute_ss_covariances(
      object1,
      fitted_values = object0$fitted.values,
      linear_predictors = object0$linear.predictors
    ),
    cov_type = cov_type,
    object = object1,
    coef_no = length(stats::coef(object1)),
    context = context
  )
}


compute_wald_test <- function(object0, object1,
                              cov_type,
                              context) {
  nested_models <- check_nested_models(object0, object1)
  obj1 <- nested_models$object1
  index <- nested_models$index
  test_df <- check_test_index(index, context)
  test_coefficients <- as.numeric(obj1$coefficients[index])
  cov_test <- stats::vcov(obj1, cov_type = cov_type)[index, index, drop = FALSE]
  test_stat <- tryCatch(
    as.numeric(
      crossprod(test_coefficients, solve(cov_test, test_coefficients))
    ),
    error = function(e) {
      stop(
        sprintf(
          "%s failed: the covariance matrix of the Wald test is singular",
          context
        ),
        call. = FALSE
      )
    }
  )
  test_stat <- check_test_statistic(test_stat, context, "Wald test")
  test_p <- stats::pchisq(test_stat, df = test_df, lower.tail = FALSE)
  list(test_stat = test_stat, test_df = test_df, test_p = test_p)
}


compute_working_wald_test <- function(object0, object1,
                                      cov_type,
                                      pmethod,
                                      context) {
  nested_models <- check_nested_models(object0, object1)
  obj1 <- nested_models$object1
  index <- nested_models$index
  check_test_index(index, context)
  test_coefficients <- as.numeric(obj1$coefficients[index])
  naive_mat <- stats::vcov(obj1, cov_type = "naive")
  cov_test <- naive_mat[index, index, drop = FALSE]
  test_stat <- tryCatch(
    as.numeric(
      crossprod(test_coefficients, solve(cov_test, test_coefficients))
    ),
    error = function(e) {
      stop(
        sprintf(
          paste0(
            "%s failed: the naive covariance matrix of the working Wald test ",
            "is singular"
          ),
          context
        ),
        call. = FALSE
      )
    }
  )
  test_stat <- check_test_statistic(test_stat, context, "working Wald test")
  robust_mat <- stats::vcov(obj1, cov_type = cov_type)
  rob_test <- robust_mat[index, index, drop = FALSE]
  eigen_test <- compute_mixture_eigenvalues(
    naive_mat = cov_test,
    robust_mat = rob_test,
    context = context,
    test_label = "working Wald test"
  )
  compute_chisq_mixture(eigen_test, test_stat, pmethod, context)
}


compute_working_lrt_test <- function(object0, object1,
                                     cov_type,
                                     pmethod,
                                     context) {
  nested_models <- check_nested_models(object0, object1)
  obj0 <- nested_models$object0
  obj1 <- nested_models$object1
  index <- nested_models$index
  check_test_index(index, context)
  phi_common <- obj1$phi
  ll_obj0 <- compute_fit_quasi_lik(obj0, phi = phi_common)
  ll_obj1 <- compute_fit_quasi_lik(obj1, phi = phi_common)
  test_stat <- -2 * (ll_obj0 - ll_obj1)
  test_stat <- check_test_statistic(test_stat, context, "working LRT")
  naive_mat <- stats::vcov(obj1, cov_type = "naive")
  cov_test <- naive_mat[index, index, drop = FALSE]
  robust_mat <- stats::vcov(obj1, cov_type = cov_type)
  rob_test <- robust_mat[index, index, drop = FALSE]
  eigen_test <- compute_mixture_eigenvalues(
    naive_mat = cov_test,
    robust_mat = rob_test,
    context = context,
    test_label = "working LRT"
  )
  compute_chisq_mixture(eigen_test, test_stat, pmethod, context)
}


compute_score_test <- function(object0, object1,
                               cov_type,
                               context) {
  nested_models <- check_nested_models(object0, object1)
  obj0 <- nested_models$object0
  obj1 <- nested_models$object1
  index <- nested_models$index
  test_df <- check_test_index(index, context)
  sc <- compute_score_components(obj0, obj1)
  score_vector <- sc$score_vector
  cov_test <- if (cov_type %in% geer_ss_cov_types) {
    compute_score_ss_cov(
      obj0, obj1, cov_type, context
    )[index, index, drop = FALSE]
  } else {
    switch(
      cov_type,
      robust = sc$robust_covariance[index, index, drop = FALSE],
      naive = sc$naive_covariance[index, index, drop = FALSE],
      `bias-corrected` = {
        check_bc_covariance(
          obj1$clusters_no,
          length(obj1$coefficients),
          context
        )
        sc$bc_covariance[index, index, drop = FALSE]
      },
      `df-adjusted` = compute_df_adjusted_covariance(
        robust_covariance = sc$robust_covariance[index, index, drop = FALSE],
        clusters_no = obj1$clusters_no,
        coef_no = length(stats::coef(obj1)),
        context = context
      ),
      jackknife = stats::vcov(obj1, cov_type = "jackknife")[
        index,
        index,
        drop = FALSE
      ]
    )
  }
  naive_mat <- sc$naive_covariance
  mid <- naive_mat[index, , drop = FALSE] %*% score_vector
  sol <- tryCatch(
    solve(cov_test, mid),
    error = function(e) {
      stop(
        sprintf(
          "%s failed: the covariance matrix of the score test is singular",
          context
        ),
        call. = FALSE
      )
    }
  )
  test_stat <- as.numeric(
    (t(score_vector) %*% naive_mat[, index, drop = FALSE]) %*% sol
  )
  test_stat <- check_test_statistic(test_stat, context, "score test")
  test_p <- stats::pchisq(test_stat, df = test_df, lower.tail = FALSE)
  list(test_stat = test_stat, test_df = test_df, test_p = test_p)
}


compute_working_score_test <- function(object0, object1,
                                       cov_type,
                                       pmethod,
                                       context) {
  nested_models <- check_nested_models(object0, object1)
  obj0 <- nested_models$object0
  obj1 <- nested_models$object1
  index <- nested_models$index
  check_test_index(index, context)
  sc <- compute_score_components(obj0, obj1)
  score_vector <- sc$score_vector
  naive_mat <- sc$naive_covariance
  robust_mat <- if (cov_type %in% geer_ss_cov_types) {
    compute_score_ss_cov(
      obj0, obj1, cov_type, context
    )
  } else {
    switch(
      cov_type,
      robust = sc$robust_covariance,
      naive = sc$naive_covariance,
      `bias-corrected` = {
        check_bc_covariance(
          obj1$clusters_no,
          length(obj1$coefficients),
          context
        )
        sc$bc_covariance
      },
      `df-adjusted` = compute_df_adjusted_covariance(
        robust_covariance = sc$robust_covariance,
        clusters_no = obj1$clusters_no,
        coef_no = length(stats::coef(obj1)),
        context = context
      ),
      jackknife = stats::vcov(obj1, cov_type = "jackknife")
    )
  }
  cov_test <- naive_mat[index, index, drop = FALSE]
  mid <- naive_mat[index, , drop = FALSE] %*% score_vector
  sol <- tryCatch(
    solve(cov_test, mid),
    error = function(e) {
      stop(
        sprintf(
          paste0(
            "%s failed: the naive covariance matrix of the working score test ",
            "is singular"
          ),
          context
        ),
        call. = FALSE
      )
    }
  )
  test_stat <- as.numeric(
    (t(score_vector) %*% naive_mat[, index, drop = FALSE]) %*% sol
  )
  test_stat <- check_test_statistic(test_stat, context, "working score test")
  eigen_test <- compute_mixture_eigenvalues(
    naive_mat = cov_test,
    robust_mat = robust_mat[index, index, drop = FALSE],
    context = context,
    test_label = "working score test"
  )
  compute_chisq_mixture(eigen_test, test_stat, pmethod, context)
}


compute_nested_test <- function(object0, object1, test, cov_type, pmethod,
                                context) {
  switch(
    test,
    wald = compute_wald_test(object0, object1, cov_type, context),
    score = compute_score_test(object0, object1, cov_type, context),
    `working-wald` = compute_working_wald_test(
      object0,
      object1,
      cov_type,
      pmethod,
      context
    ),
    `working-score` = compute_working_score_test(
      object0,
      object1,
      cov_type,
      pmethod,
      context
    ),
    `working-lrt` = compute_working_lrt_test(
      object0,
      object1,
      cov_type,
      pmethod,
      context
    )
  )
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
    stop(
      sprintf(
        paste0(
          "models must have the same response; differing response(s): %s"
        ),
        paste(removed, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  association_vector <- vapply(
    object,
    function(x) x$association_structure,
    character(1)
  )
  for (structure_name in association_vector) {
    check_working_lrt_allowed(test, structure_name)
  }
  same_association <- association_vector == association_vector[1L]
  if (!all(same_association)) {
    removed <- unique(association_vector[!same_association])
    stop(
      sprintf(
        paste0(
          "models must have the same working association structure; ",
          "differing structure(s): %s"
        ),
        paste(removed, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  sample_size_vector <- vapply(object, get_model_sample_size, numeric(1))
  if (any(sample_size_vector != sample_size_vector[1L])) {
    stop(
      "models must be fitted to the same observations",
      call. = FALSE
    )
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
  test_type <- build_test_label(test)
  for (i in 2:models_no) {
    value <- compute_nested_test(
      object[[i - 1L]],
      object[[i]],
      test,
      cov_type,
      pmethod,
      "anova.geer"
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
