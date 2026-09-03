validate_mcar_little_matrix <- function(x) {
  if (is.data.frame(x)) {
    numeric_cols <- vapply(x, is.numeric, logical(1))
    if (!all(numeric_cols)) {
      stop("all variables supplied to Little's MCAR test must be numeric", call. = FALSE)
    }
    x <- as.matrix(x)
  } else if (is.matrix(x)) {
    if (!is.numeric(x)) {
      stop("the matrix supplied to Little's MCAR test must be numeric", call. = FALSE)
    }
  } else {
    stop("'object' must be a fitted 'geer' object, numeric matrix, or numeric data frame", call. = FALSE)
  }

  storage.mode(x) <- "double"
  if (nrow(x) < 2L) {
    stop("Little's MCAR test requires at least two rows", call. = FALSE)
  }
  if (ncol(x) < 2L) {
    stop("Little's MCAR test requires at least two variables or repeated measurements", call. = FALSE)
  }
  if (any(is.infinite(x))) {
    stop("data supplied to Little's MCAR test cannot contain infinite values", call. = FALSE)
  }
  observed_no <- colSums(!is.na(x))
  labels <- colnames(x)
  if (is.null(labels)) labels <- as.character(seq_len(ncol(x)))
  if (any(observed_no < 2L)) {
    bad <- which(observed_no < 2L)
    stop(
      sprintf(
        "Little's MCAR test requires at least two observed values for variable(s) %s",
        paste(labels[bad], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (ncol(x) > 1L) {
    pairs <- utils::combn(ncol(x), 2L)
    overlap <- vapply(
      seq_len(ncol(pairs)),
      function(j) {
        sum(!is.na(x[, pairs[1L, j]]) & !is.na(x[, pairs[2L, j]]))
      },
      integer(1)
    )
    if (any(overlap < 1L)) {
      bad_pair <- pairs[, which(overlap < 1L)[1L]]
      stop(
        sprintf(
          paste0(
            "Little's MCAR test requires at least one jointly observed value ",
            "for every variable pair; '%s' and '%s' do not meet this requirement"
          ),
          labels[bad_pair[1L]],
          labels[bad_pair[2L]]
        ),
        call. = FALSE
      )
    }
  }
  x
}


mcar_little_bivariate_exact <- function(x) {
  n <- nrow(x)
  if (ncol(x) != 2L || n <= 2L) return(NULL)

  observed_no <- colSums(!is.na(x))
  complete_col <- which(observed_no == n)
  incomplete_col <- which(observed_no < n)
  if (length(complete_col) != 1L || length(incomplete_col) != 1L) return(NULL)

  missing_group <- is.na(x[, incomplete_col])
  if (!any(missing_group) || all(missing_group)) return(NULL)

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

  # Little (1988), Section 3.4: d^2 = (n - 1) F / (n - 2 + F).
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


mcar_little_calculate <- function(x, maxit, tol) {
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

  em <- mcar_normal_em(x, maxit = maxit, tol = tol)
  sigma_corrected <- (n / (n - 1)) * em$sigma
  statistic <- 0
  df_terms <- 0L

  for (rows in pattern_rows) {
    observed <- which(!missing[rows[1L], ])
    k <- length(observed)
    df_terms <- df_terms + k
    if (!k) next

    pattern_mean <- colMeans(x[rows, observed, drop = FALSE])
    difference <- pattern_mean - em$mu[observed]
    sigma_observed <- sigma_corrected[observed, observed, drop = FALSE]
    solved <- tryCatch(
      solve(sigma_observed, difference),
      error = function(e) {
        stop(
          "a covariance submatrix required for Little's MCAR statistic is singular",
          call. = FALSE
        )
      }
    )
    statistic <- statistic + length(rows) * as.numeric(crossprod(difference, solved))
  }

  df <- as.integer(df_terms - p)
  if (df <= 0L) {
    stop(
      "Little's MCAR test has non-positive degrees of freedom for this missing-data pattern",
      call. = FALSE
    )
  }

  statistic <- max(as.numeric(statistic), 0)
  asymptotic_p_value <- stats::pchisq(statistic, df = df, lower.tail = FALSE)
  bivariate_exact <- mcar_little_bivariate_exact(x)

  if (!is.null(bivariate_exact) &&
      !isTRUE(all.equal(
        statistic,
        bivariate_exact$implied_little_statistic,
        tolerance = 1e-6,
        check.attributes = FALSE
      ))) {
    warning(
      paste0(
        "the bivariate monotone small-sample identity in Little (1988) was ",
        "not reproduced to numerical tolerance; using the asymptotic reference"
      ),
      call. = FALSE
    )
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


mcar_little_has_binary_variable <- function(x) {
  any(vapply(
    seq_len(ncol(x)),
    function(j) {
      values <- unique(x[!is.na(x[, j]), j])
      length(values) >= 2L && all(values %in% c(0, 1))
    },
    logical(1)
  ))
}

