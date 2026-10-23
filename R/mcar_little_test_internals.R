check_mcar_little_matrix <- function(x) {
  x <- check_mcar_matrix(x, "Little's MCAR test")
  labels <- colnames(x)
  if (is.null(labels)) {
    labels <- as.character(seq_len(ncol(x)))
  }
  pairs <- utils::combn(ncol(x), 2L)
  overlap <- vapply(
    seq_len(ncol(pairs)),
    function(j) {
      sum(!is.na(x[, pairs[1L, j]]) & !is.na(x[, pairs[2L, j]]))
    },
    numeric(1)
  )
  if (any(overlap < 1L)) {
    bad_pair <- pairs[, which(overlap < 1L)[1L]]
    stop(
      sprintf(
        paste0(
          "Little's MCAR test requires every pair of variables to be ",
          "observed together at least once; '%s' and '%s' never are"
        ),
        labels[bad_pair[1L]],
        labels[bad_pair[2L]]
      ),
      call. = FALSE
    )
  }
  x
}


compute_little_exact <- function(x) {
  n <- nrow(x)
  if (ncol(x) != 2L || n <= 2L) {
    return(NULL)
  }

  observed_no <- colSums(!is.na(x))
  complete_col <- which(observed_no == n)
  incomplete_col <- which(observed_no < n)
  if (length(complete_col) != 1L || length(incomplete_col) != 1L) {
    return(NULL)
  }

  missing_group <- is.na(x[, incomplete_col])
  if (!any(missing_group) || all(missing_group)) {
    return(NULL)
  }

  y <- x[, complete_col]
  group_mean_observed <- mean(y[!missing_group])
  group_mean_missing <- mean(y[missing_group])
  overall_mean <- mean(y)
  n_observed <- sum(!missing_group)
  n_missing <- sum(missing_group)

  ss_between <- n_observed * (group_mean_observed - overall_mean)^2 +
    n_missing * (group_mean_missing - overall_mean)^2
  ss_within <- sum((y[!missing_group] - group_mean_observed)^2) +
    sum((y[missing_group] - group_mean_missing)^2)
  df2 <- n - 2L
  f_statistic <- if (ss_within == 0) Inf else ss_between / (ss_within / df2)
  p_value <- stats::pf(f_statistic, df1 = 1, df2 = df2, lower.tail = FALSE)

  implied_statistic <- if (is.infinite(f_statistic)) {
    n - 1
  } else {
    (n - 1) * f_statistic / (n - 2 + f_statistic)
  }

  list(
    statistic = f_statistic,
    df1 = 1L,
    df2 = as.integer(df2),
    p_value = p_value,
    implied_little_statistic = implied_statistic
  )
}


compute_mcar_little_statistic <- function(x, maxiter, tolerance, reference) {
  missing <- is.na(x)
  n <- nrow(x)
  p <- ncol(x)
  pattern_key <- apply(missing, 1L, paste0, collapse = "")
  pattern_rows <- split(seq_len(n), pattern_key)
  missing_patterns <- length(pattern_rows)

  if (!any(missing)) {
    sigma_ml <- stats::cov(x) * (n - 1) / n
    sigma_corrected <- stats::cov(x)
    return(list(
      statistic = 0,
      df = 0L,
      p_value = 1,
      missing_patterns = 1L,
      iterations = 0L,
      converged = TRUE,
      mu = colMeans(x),
      sigma = sigma_corrected,
      sigma_ml = sigma_ml,
      bivariate_exact = NULL
    ))
  }

  em <- fit_mcar_normal_em(
    x,
    maxiter = maxiter,
    tolerance = tolerance,
    context = "mcar_little_test"
  )
  normalized_x <- em$normalized$values
  normalized_sigma <- (n / (n - 1)) * em$normalized_sigma
  sigma_corrected <- (n / (n - 1)) * em$sigma
  statistic <- 0
  df_terms <- 0L

  for (rows in pattern_rows) {
    observed <- which(!missing[rows[1L], ])
    df_terms <- df_terms + length(observed)

    pattern_mean <- colMeans(normalized_x[rows, observed, drop = FALSE])
    difference <- pattern_mean - em$normalized_mu[observed]
    sigma_observed <- normalized_sigma[observed, observed, drop = FALSE]
    chol_factor <- tryCatch(chol(sigma_observed), error = function(e) NULL)
    if (is.null(chol_factor)) {
      stop(
        paste0(
          "mcar_little_test failed: a covariance submatrix required for ",
          "Little's MCAR statistic is not positive definite"
        ),
        call. = FALSE
      )
    }
    solved <- forwardsolve(t(chol_factor), difference)
    statistic <- statistic + length(rows) * sum(solved^2)
  }

  df <- as.integer(df_terms - p)
  if (df <= 0L) {
    stop(
      paste0(
        "mcar_little_test failed: Little's MCAR test has non-positive ",
        "degrees of freedom for this missing-data pattern"
      ),
      call. = FALSE
    )
  }

  statistic <- as.numeric(statistic)
  asymptotic_p_value <- stats::pchisq(statistic, df = df, lower.tail = FALSE)
  bivariate_exact <- compute_little_exact(normalized_x)

  if (!is.null(bivariate_exact) &&
      !isTRUE(all.equal(
        statistic,
        bivariate_exact$implied_little_statistic,
        tolerance = 1e-6,
        check.attributes = FALSE
      ))) {
    if (identical(reference, "auto")) {
      stop(
        paste0(
          "mcar_little_test failed: the small-sample identity of Little ",
          "(1988) is not reproduced numerically; use reference = ",
          "\"asymptotic\""
        ),
        call. = FALSE
      )
    }
    bivariate_exact <- NULL
  }

  list(
    statistic = statistic,
    df = df,
    p_value = asymptotic_p_value,
    missing_patterns = as.integer(missing_patterns),
    iterations = as.integer(em$iterations),
    converged = isTRUE(em$converged),
    mu = em$mu,
    sigma = sigma_corrected,
    sigma_ml = em$sigma,
    bivariate_exact = bivariate_exact
  )
}


is_mcar_little_binary <- function(x) {
  any(vapply(
    seq_len(ncol(x)),
    function(j) {
      values <- unique(x[!is.na(x[, j]), j])
      length(values) >= 2L && all(values %in% c(0, 1))
    },
    logical(1)
  ))
}
